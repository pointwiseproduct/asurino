#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <memory>
#include <string>
#include <sstream>
#include <functional>
#include <regex>
#include <random>
#include <chrono>
#include <thread>
#include <cmath>
#include <cstdlib>
#include <mecab.h>

std::unique_ptr<MeCab::Model> mecab_model(MeCab::createModel(""));
std::unique_ptr<MeCab::Tagger> mecab_tagger(MeCab::createTagger(""));

double string_to_double(const std::string &str){
    double r = 0.0;
    int phase = 0;
    for(std::size_t i = 0; i < str.size(); ++i){
        if(str[i] == '.'){
            phase = 1;
            continue;
        }
        int n = str[i] - '0';
        if(phase == 0){
            r *= 10;
            r += n;
        }else{
            r += std::pow(10.0, -phase) * n;
            ++phase;
        }
    }
    return r;
}

struct rule{
    rule() = default;
    rule(const rule&) = default;
    rule(rule&&) = default;
    ~rule() = default;

    rule(const std::string &expression){
        *this = form_string(expression);
    }

    rule &operator =(const rule &rhs){
        source = rhs.source;
        result1 = rhs.result1;
        result2 = rhs.result2;
        probability = rhs.probability;
        return *this;
    }

    bool operator ==(const rule &other) const{
        return source == other.source && result1 == other.result1 && result2 == other.result2;
    }

    std::string source;
    std::string result1;
    std::string result2;
    double probability;

    static rule form_string(const std::string &expression){
        rule r;
        std::regex pattern(R"((.+)\s*>\s*(.+)\s+(.+)\s+(.+))");
        std::cmatch result;
        if(std::regex_match(expression.c_str(), result, pattern)){
            r.source = result.str(1);
            r.result1 = result.str(2);
            r.result2 = result.str(3);
            r.probability = string_to_double(result.str(4));
        }
        return r;
    }

    static rule end_rule(const std::string &pos){
        rule r;
        r.source = pos;
        r.result1 = "END";
        r.result2 = "END";
        r.probability = 1.0;
        return r;
    }

    std::string to_string() const{
        return source + ">" + result1 + " " + result2 + " " + std::to_string(probability);
    }
};

namespace rule_tree{
    struct coordinate{
        double x;
        double y;

        bool operator ==(const coordinate &other) const{
            return x == other.x && y == other.y;
        }
    };

    using rule_coordinate_map = std::map<std::string, coordinate>;
    using pos_rule_list_map = std::map<std::string, std::vector<rule>>;

    struct node{
        rule_coordinate_map left, right;
        pos_rule_list_map rules;

        bool operator ==(const node &other) const{
            return left == other.left && right == other.right && rules == other.rules;
        }
    };
}

struct token;
using token_vector = std::vector<token>;

struct token{
    // +1分はrawな形態素データ
    static const std::size_t size = 9 + 1;

    token(){
        for(int i = 0; i < size; ++i){
            data[i] = "*";
        }
    }

    token(const token&) = default;

    token(token &&other){
        for(int i = 0; i < sizeof(data) / sizeof(data[0]); ++i){
            data[i] = std::move(other.data[i]);
        }
    }

    ~token() = default;

    token(const MeCab::Node &n){
        input(n);
    }

    token &operator =(const token &other){
        for(int i = 0; i < size; ++i){
            data[i] = other.data[i];
        }
        return *this;
    }
    
    std::string data[size];
    const std::string *data_end = data + size;

    // 品詞．
    std::string &noun = data[0];

    // 品詞の詳細．
    std::string *noun_detail = &data[0];
    std::string *noun_detail_end = &data[6];

    // 形態素の形．
    std::string &form = data[5];

    // 形態素の文字列．
    std::string &str = data[size - 1];

    bool operator <(const token &other) const{
        for(int i = 0; i < size; ++i){
            if(data[i] < other.data[i]){
                return true;
            }else if(data[i] > other.data[i]){
                return false;
            }
        }
        return false;
    }

    bool operator ==(const token &other) const{
        for(int i = 0; i < size; ++i){
            if(data[i] != other.data[i]){
                return false;
            }
        }
        return true;
    }

    void clear_unique_data(){
        for(std::size_t i = 6; i < size; ++i){
            data[i] = "*";
        }
    }

    bool is_eos() const{
        static const std::string str = "BOS/EOS";
        return noun == str;
    }

    void input(const MeCab::Node &n){
        const char *raw_data = n.feature;

        for(int i = 0; i < size; ++i){
            data[i].clear();
        }

        data[size - 1] = std::string(n.surface, n.surface + n.length);

        std::size_t str_count = 0;
        for(int count = 0; raw_data[str_count]; ++count){
            while(true){
                if(raw_data[str_count] != ',' && raw_data[str_count]){
                    data[count] += raw_data[str_count];
                    ++str_count;
                }else{
                    if(raw_data[str_count] == ','){
                        ++str_count;
                    }
                    break;
                }
            }
        }
    }

    // 品詞情報に該当するものがあるかどうか調べる
    bool contain_noun_detail(const std::string &info) const{
        for(std::size_t i = 1; i <= 5; ++i){
            if(data[i] == info){
                return true;
            }
        }
        return false;
    }

    // 文字列にシリアライズする
    std::string serialize() const{
        std::string r;
        for(std::size_t i = 0; i < size; ++i){
            if(i){
                r += ",";
            }
            r += data[i];
        }
        return r;
    }

    // シリアライズされた文字列から復元する
    void restore(const std::string str){
        std::size_t count = 0;
        for(std::size_t i = 0; i < size; ++i){
            data[i] = "";
            data[i] += str[count];
            ++count;
            while(count < str.size() && str[count] != ','){
                data[i] += str[count];
                ++count;
            }
            if(count >= str.size()){
                break;
            }
        }
    }

    // 品詞マッチ．
    // 品詞情報 "*" は無制限にマッチする．
    static bool noun_match(const token &lhs, const token &rhs){
        std::string *liter = lhs.noun_detail;
        std::string *riter = lhs.noun_detail;
        for(; liter != lhs.noun_detail_end; ++liter, ++riter){
            if(*liter == "*" || *riter == "*"){
                continue;
            }else if(*liter == *riter){
                continue;
            }else{
                return false;
            }
        }
        return true;
    }

    static token_vector make_vector(const MeCab::Node *first){
        token_vector r;
        const MeCab::Node *iter = first->next;
        while(true){
            token t(*iter);
            if(t.is_eos()){
                break;
            }
            r.push_back(std::move(t));
            iter = iter->next;
        }
        return r;
    }
};

class pcfg{
public:
    using pos_map = std::map<std::string, double>;

    class pcfg_node{
    public:
        pos_map inside, outside;
    };

    static void parse(
        const std::string &text,
        const std::vector<rule> rules,
        std::function<void(
            const std::vector<std::vector<rule_tree::node>>*,
            const std::vector<token>&,
            const std::vector<rule>&
        )> callback
    ){
        pcfg parser;
        auto fn = [&](const std::vector<token> &tokens){
            parser.rules = rules;
            bool parsed = parser.calc(tokens);
            if(parsed){
                parser.recalc_probability(tokens, parser.node_map);
                callback(&parser.rule_tree_map, tokens, parser.rules);
            }else{
                callback(nullptr, tokens, parser.rules);
            }
        };
        parser.tokenize(text, fn);
    }

    void tokenize(const std::string &text, std::function<void(const std::vector<token>&)> callback){
        std::vector<token> tokens = token::make_vector(mecab_tagger->parseToNode(text.c_str()));
        callback(tokens);
    }

    bool calc(const std::vector<token> &tokens){
        std::size_t N = tokens.size();
        node_map.resize(N);
        rule_tree_map.resize(N);
        for(std::size_t i = 0; i < N; ++i){
            node_map[i].resize(N);
            rule_tree_map[i].resize(N);
        }
        for(std::size_t i = 0; i < N; ++i){
            const std::string &pos = tokens[i].noun;
            node_map[i][i].inside[pos] = 1.0;
            rule_tree_map[i][i].rules[pos].resize(1);
            rule_tree_map[i][i].rules[pos][0] = rule::end_rule(pos);
        }
        // ボトムアップ．
        for(std::size_t n = 1; n < N; ++n){
            for(std::size_t i = 0; i < N - n; ++i){
                std::size_t x = i;
                std::size_t y = i + n;
                for(std::size_t j = 1; j < 1 + n; ++j){
                    for(std::size_t rule_index = 0; rule_index < rules.size(); ++rule_index){
                        rule &rule_ = rules[rule_index];
                        // CYK．
                        double *left_ptr = nullptr;
                        {
                            auto &map = node_map[x][x + (j - 1)].inside;
                            auto iter = map.find(rule_.result1);
                            if(iter == map.end()){
                                iter = map.insert(iter, std::make_pair(rule_.result1, 0));
                            }
                            left_ptr = &iter->second;
                        }
                        const double &left = *left_ptr;

                        double *bottom_ptr = nullptr;
                        {
                            auto &map = node_map[x + j][y].inside;
                            auto iter = map.find(rule_.result2);
                            if(iter == map.end()){
                                iter = map.insert(iter, std::make_pair(rule_.result2, 0));
                            }
                            bottom_ptr = &iter->second;
                        }
                        const double &bottom = *bottom_ptr;

                        double *old_ptr = nullptr;
                        {
                            auto &map = node_map[x][y].inside;
                            auto iter = map.find(rule_.source);
                            if(iter == map.end()){
                                iter = map.insert(iter, std::make_pair(rule_.source, 0));
                            }
                            old_ptr = &iter->second;
                        }
                        const double &old = *old_ptr;

                        node_map[x][y].inside[rule_.source] = old + rule_.probability * left * bottom;

                        if(rule_.probability * left * bottom > 0){
                            // 最尤推定用．
                            rule_tree_map[x][y].rules[rule_.source].push_back(rule_);
                            rule_tree::coordinate &left = rule_tree_map[x][y].left[rule_.to_string()];
                            left.x = x;
                            left.y = x + (j - 1);

                            rule_tree::coordinate &right = rule_tree_map[x][y].right[rule_.to_string()];
                            right.x = x + j;
                            right.y = y;
                        }
                    }
                }
            }
        }
        // 非受理判定．
        if(node_map[0][N - 1].inside["S"] == 0){
            return false;
        }
        // 外側初期化．
        for(std::size_t i = 0; i < N; ++i){
            for(std::size_t n = 0; n < N - i; ++n){
                for(auto &source : to_unique_array(map<rule, std::string>(rules, [](const rule &f){ return f.source; }))){
                    if(node_map[i][i + n].inside[source] == 0){
                        node_map[i][i + n].outside[source] = 0;
                    }
                }
            }
        }
        // 頂点定義．
        node_map[0][N - 1].outside["S"] = 1;
        // トップダウン．
        for(int n = N - 1; n >= 0; --n){
            for(std::size_t i = 0; i < N - n; ++i){
                std::size_t x = i;
                std::size_t y = i + n;
                for(std::size_t rule_index = 0; rule_index < rules.size(); ++rule_index){
                    const rule &rule_ = rules[rule_index];
                    if(node_map[x][y].outside[rule_.source] > 0){
                        for(std::size_t k = 1; static_cast<int>(k) < 1 + n; ++k){
                            // 下側．
                            if(node_map[x + k][y].inside[rule_.result2] > 0){
                                double value = 
                                    node_map[x + k][y].inside[rule_.result2] +
                                    (rule_.probability *
                                        node_map[x][y].outside[rule_.source] *
                                        node_map[x][x - (1 - k)].inside[rule_.result1]);
                                node_map[x + k][y].outside[rule_.result2] = value;
                            }
                            // 左側．
                            if(node_map[x][y - k].inside[rule_.result1] > 0){
                                double value =
                                    node_map[x][y - k].outside[rule_.result1] +
                                    (rule_.probability *
                                        node_map[x][y].outside[rule_.source] *
                                        node_map[y + (1 - k)][y].inside[rule_.result2]);
                                node_map[x][y - k].outside[rule_.result1] = value;
                            }
                        }
                    }
                }
            }
        }
        return true;
    }

    void recalc_probability(const std::vector<token> &tokens, std::vector<std::vector<pcfg_node>> &node_map){
        auto filter = [](const std::vector<rule> &rules, std::function<bool(const rule&)> f){
            std::vector<rule> r;
            for(auto &i : rules){
                if(f(i)){
                    r.push_back(i);
                }
            }
            return r;
        };

        auto reduce = [](const std::vector<rule> &rules, std::function<double(const double&, const rule&)> f, double r){
            if(rules.size() == 1){
                return f(r, rules[0]);
            }
            for(std::size_t i = 1; i < rules.size(); ++i){
                r = f(r, rules[i]);
            }
            return r;
        };

        for(std::size_t i = 0; i < rules.size(); ++i){
            rule &rule_ = rules[i];
            double new_p
                = use_count(tokens.size(), rule_, node_map)
                / reduce(filter(rules, [&](const rule &x){ return x.source == rule_.source; }), [&](const double &pv, const rule &x){ return pv + use_count(tokens.size(), x, node_map); }, 0);
            rule_.probability = (rule_.probability + new_p) / 2;
        }
    }

    static double use_count(std::size_t N, const rule &rule_, std::vector<std::vector<pcfg_node>> &node_map){
        std::size_t count = 0;
        for(std::size_t n = 1; n < N; ++n){
            for(std::size_t i = 0; i < N - n; ++i){
                std::size_t x = i;
                std::size_t y = i + n;
                for(std::size_t j = 1; j < n + 1; ++j){
                    double source;
                    {
                        auto &map = node_map[x][y].outside;
                        auto iter = map.find(rule_.source);
                        if(iter == map.end()){
                            iter = map.insert(std::make_pair(rule_.source, 0)).first;
                        }
                        source = iter->second;
                    }

                    double result1;
                    {
                        auto &map = node_map[x][x + j - 1].inside;
                        auto iter = map.find(rule_.result1);
                        if(iter == map.end()){
                            iter = map.insert(std::make_pair(rule_.result1, 0)).first;
                        }
                        result1 = iter->second;
                    }

                    double result2;
                    {
                        auto &map = node_map[x + j][y].inside;
                        auto iter = map.find(rule_.result2);
                        if(iter == map.end()){
                            iter = map.insert(std::make_pair(rule_.result2, 0)).first;
                        }
                        result2 = iter->second;
                    }

                    count += source * result1 * result2 > 0 ? 1 : 0;
                }
            }
        }
        return count <= 0 ? 000000000000001 : count;
    }

private:
    template<class T>
    std::vector<T> to_unique_array(const std::vector<T> &array){
        auto indexof = [](const std::vector<T> &arr, const T &t) -> int{
            for(int i = 0; i < static_cast<int>(arr.size()); ++i){
                if(arr[i] == t){
                    return i;
                }
            }
            return -1;
        };

        std::vector<T> a;
        for(std::size_t i = 0, l = array.size(); i < l; ++i){
            if(indexof(a, array[i]) == -1){
                a.push_back(array[i]);
            }
        }
        return a;
    }

    template<class T, class U>
    std::vector<U> map(const std::vector<T> &v, U f(const T&)){
        std::vector<U> r;
        r.reserve(v.size());
        for(auto &i : v){
            r.push_back(f(i));
        }
        return r;
    }

    std::vector<rule> rules;
    std::vector<std::vector<pcfg_node>> node_map;
    std::vector<std::vector<rule_tree::node>> rule_tree_map;
};

void test_myparser(){
    std::vector<rule> rules = {
        rule("S>名詞句 形容詞 0.5"),
        rule("S>名詞句 名詞 0.3"),
        rule("S>名詞句 動詞 0.2"),
        rule("名詞>形容詞 名詞 1.0"),
        rule("名詞句>名詞 助詞 1.0"),
        rule("形容詞>名詞 助詞 0.1"),
        rule("形容詞>名詞 動詞 0.2"),
        rule("形容詞>副詞 形容詞 0.4"),
        rule("形容詞>副詞 形容詞 0.3"),
        rule("動詞>副詞 動詞句 0.5"),
        rule("動詞>名詞 助動詞 0.5")
    };

    std::function<std::size_t(
        const std::vector<std::vector<rule_tree::node>>&,
        const std::vector<token>&,
        std::size_t x, std::size_t,
        const std::string&,
        std::size_t,
        std::size_t
    )> display;
    display = [&](
        const std::vector<std::vector<rule_tree::node>> &tree,
        const std::vector<token> &tokens,
        std::size_t x, std::size_t y,
        const std::string &pos,
        std::size_t depth, std::size_t leaf_count
    ){
        const rule_tree::node &top = tree[x][y];
        if(top == rule_tree::node()){
            return leaf_count;
        }
        std::string result = "";
        std::vector<rule> temp = top.rules.find(pos)->second;
        std::sort(temp.begin(), temp.end(), [](const rule &l, const rule &r){ return l.probability < r.probability; });
        rule rule_ = temp[0];
        if(rule_.result1 == "END"){
            result = "---> " + tokens[leaf_count].str;
            ++leaf_count;
        }else{
            result = "(" + std::to_string(rule_.probability) + ")";
        }

        for(std::size_t i = 0; i < depth; ++i){
            std::cout << "    ";
        }
        std::cout << rule_.source << result << std::endl;

        if(rule_.result1 != "END"){
            leaf_count = display(
                tree,
                tokens,
                top.left.find(rule_.to_string())->second.x,
                top.left.find(rule_.to_string())->second.y,
                rule_.result1,
                depth + 1,
                leaf_count
            );
        }
        if(rule_.result2 != "END"){
            leaf_count = display(
                tree,
                tokens,
                top.right.find(rule_.to_string())->second.x,
                top.right.find(rule_.to_string())->second.y,
                rule_.result2,
                depth + 1,
                leaf_count
            );
        }
        return leaf_count;
    };
}

class word_graph{
public:
    static const std::size_t words_upper_limit = 2048;
    static const std::size_t link_upper_limit = 192;

    struct node_frame{
        token t;
        std::list<token> link;
        std::map<token, std::list<token>::iterator> link_map;

        // リンクにワードを追加する
        void add_link(const token &t){
            auto iter = link_map.find(t);
            if(iter == link_map.end()){
                link.push_front(t);
                link_map.insert(std::make_pair(t, link.begin()));
                if(link.size() >= link_upper_limit){
                    link_map.erase(link.back());
                    link.pop_back();
                }
            }
        }

        // シリアライズ
        std::string serialize() const{
            std::string r = t.serialize();
            for(auto &i : link){
                r += "\n";
                r += i.serialize();
            }
            return r;
        }

        // 複合する
        void restore(const std::string &str){
            std::stringstream ss(str);
            std::string line;
            std::getline(ss, line);
            t.restore(line);
            while(true){
                std::getline(ss, line);
                if(ss.fail()){
                    break;
                }
                token u;
                u.restore(line);
                link.push_back(u);
                auto iter = link.end();
                --iter;
                link_map.insert(std::make_pair(u, iter));
            }
        }
    };

    std::list<node_frame> data_vec;
    std::map<token, std::list<node_frame>::iterator> data_map;

    // 文字列を解析する
    void set_str(const token_vector &vec){
        std::vector<decltype(data_map)::iterator> iter_vec;
        std::vector<bool> flag_vec;
        flag_vec.reserve(vec.size());
        for(auto &i : vec){
            if(unique_word_filter(i)){
                iter_vec.push_back(repush(i));
                flag_vec.push_back(true);
            }else{
                iter_vec.push_back(data_map.end());
                flag_vec.push_back(false);
            }
        }
        for(std::size_t i = 0; i < vec.size(); ++i){
            if(!flag_vec[i]){
                continue;
            }
            for(std::size_t j = 0; j < vec.size(); ++j){
                if(!flag_vec[j]){
                    continue;
                }
                if(i == j){
                    continue;
                }
                if(iter_vec[i]->second->t == iter_vec[j]->first){
                    continue;
                }
                auto &frame = *iter_vec[i]->second;
                frame.add_link(iter_vec[j]->first);
            }
        }
    }

    // 関連語を取得する
    void get_relational_words(std::set<token> &set, const token &t, int chain) const{
        if(chain <= 0){
            return;
        }
        set.insert(t);
        for(auto &i : data_map){
            for(auto &j : i.second->link){
                set.insert(j);
                get_relational_words(set, j, chain - 1);
            }
        }
    }

    // ユニークワードをフィルターする
    static bool unique_word_filter(const token &t){
        if(t.data[6] == ""){
            return false;
        }
        if(t.data[0] == "名詞"){
            if(t.data[2] == "一般"){
                return false;
            }
            if(t.contain_noun_detail("副詞可能")){
                return false;
            }
            if(t.contain_noun_detail("助数詞")){
                return false;
            }
            if(t.contain_noun_detail("数")){
                return false;
            }
            return true;
        }
        if(t.data[0] == "動詞"){
            if(t.contain_noun_detail("非自立")){
                return false;
            }
            if(t.contain_noun_detail("一段") && t.contain_noun_detail("未然形")){
                return false;
            }
            return true;
        }
        return false;
    }

    // ワードを記録する
    // 記録内に既にあれば削除待ちキューの前にもっていく
    // (前にあると削除されない)
    std::map<token, std::list<node_frame>::iterator>::iterator repush(const token &t){
        auto iter = data_map.find(t);
        if(iter == data_map.end()){
            node_frame f;
            f.t = t;
            data_vec.push_front(f);
            iter = data_map.insert(iter, std::make_pair(t, data_vec.begin()));
            if(data_vec.size() >= words_upper_limit){
                const node_frame &g = data_vec.back();
                data_map.erase(g.t);
                data_vec.pop_back();
            }
        }else{
            std::swap(data_vec.front(), *iter->second);
            iter->second = data_vec.begin();
        }
        return iter;
    }

    void serialize(std::ostream &os) const{
        for(auto &i : data_vec){
            os << i.serialize() + "\n----\n";
        }
        os << "EOS\n";
    }

    void restore(std::istream &is){
        while(true){
            std::string line;
            std::string chunk;
            while(true){
                std::getline(is, line);
                if(is.fail() || line == "EOS"){
                    break;
                }
                if(line == "----"){
                    break;
                }
                chunk += line + "\n";
            }
            if(is.fail() || line == "EOS"){
                break;
            }
            node_frame f;
            f.restore(chunk);
            data_vec.push_back(f);
            auto iter = data_vec.end();
            --iter;
            data_map.insert(std::make_pair(f.t, iter));
        }
    }
};

class sentences{
public:
    static const std::size_t log_upper_limit = 1024;

    struct frame{
        struct node{
            token t;
            bool knockout;
        };

        std::vector<node> vec;
    };

    void add_sentence(const token_vector &raw_sentence){
        frame f;
        f.vec.reserve(raw_sentence.size());
        for(std::size_t i = 0; i < raw_sentence.size(); ++i){
            frame::node n;
            n.t = raw_sentence[i];
            if(word_graph::unique_word_filter(raw_sentence[i])){
                n.knockout = true;
            }else{
                n.knockout = false;
            }
            f.vec.push_back(n);
        }
        sentence_vec.push_back(f);
        if(sentence_vec.size() >= log_upper_limit){
            sentence_vec.erase(sentence_vec.begin(), sentence_vec.begin() + (sentence_vec.size() - log_upper_limit) + 1);
        }
    }

    void serialize(std::ostream &os) const{
        for(auto &f : sentence_vec){
            for(auto &n : f.vec){
                if(n.knockout){
                    os << "- ";
                }else{
                    os << "+ ";
                }
                os << n.t.serialize() << "\n";
            }
            os << "----\n";
        }
        os << "EOS\n";
    }

    void restore(std::istream &is){
        while(true){
            std::string line;
            std::vector<std::string> chunk;
            while(true){
                std::getline(is, line);
                if(is.fail() || line == "EOS"){
                    break;
                }
                if(line == "----"){
                    break;
                }
                chunk.push_back(line);
            }
            if(is.fail() || line == "EOS"){
                break;
            }
            frame f;
            for(const auto &s : chunk){
                frame::node n;
                if(s[0] == '-'){
                    n.knockout = true;
                }else if(s[1] == '+'){
                    n.knockout = false;
                }
                n.t.restore(s.substr(2, s.size() - 2));
                f.vec.push_back(n);
            }
            if(!f.vec.empty()){
                sentence_vec.push_back(f);
            }
        }
    }

    std::vector<frame> sentence_vec;
};

std::string user_name;

std::vector<std::string> read_log(std::istream &is){
    std::vector<std::string> raw;
    while(!is.fail()){
        std::string line;
        std::getline(is, line);
        raw.push_back(line);
    }
    const std::string head = "[0m[92;49m";
    std::vector<std::string> content;
    bool asurinos_post = false;
    for(auto &i : raw){
        if(i.size() > head.size() && i.substr(0, head.size()) == head){
            asurinos_post = i.find(user_name) != std::string::npos;
            continue;
        }
        if(!asurinos_post){
            content.push_back(i);
        }
    }

    std::string head2;
    head2.push_back(27);
    head2.push_back(91);
    head2.push_back(48);
    head2.push_back(109);
    head2.push_back(27);
    head2.push_back(91);
    head2.push_back(51);
    head2.push_back(55);
    head2.push_back(59);
    head2.push_back(52);
    head2.push_back(57);
    head2.push_back(109);

    std::string foot;
    foot.push_back(27);
    foot.push_back(91);
    foot.push_back(48);
    foot.push_back(109);

    std::vector<std::string> content2;
    std::size_t count = 0;
    std::string tw;
    content[0] = head2;
    while(count < content.size()){
        std::string str = content[count];
        if(str.size() >= head2.size() && str.substr(0, head2.size()) == head2){
            if(!tw.empty()){
                std::string str1 = std::regex_replace(tw, std::regex(R"(http(s)?://([\w-]+\.)+[\w-]+(/[\w- ./?%&=]*)?)"), "");
                std::string str2 = std::regex_replace(str1, std::regex(R"(@[a-zA-Z0-9_]+)"), "");
                content2.push_back(str2);
            }
            if(str.size() > foot.size() && str.substr(str.size() - foot.size(), foot.size()) == foot){
                str = str.substr(0, str.size() - foot.size());
            }
            tw = str.substr(head2.size(), str.size() - head2.size());
        }else{
            if(str.size() > foot.size() && str.substr(str.size() - foot.size(), foot.size()) == foot){
                str = str.substr(0, str.size() - foot.size());
            }
            tw += "\n" + str;
        }
        ++count;
    }
    if(!tw.empty()){
        std::string str1 = std::regex_replace(tw, std::regex(R"(http(s)?://([\w-]+\.)+[\w-]+(/[\w- ./?%&=]*)?)"), "");
        std::string str2 = std::regex_replace(str1, std::regex(R"(@[a-zA-Z0-9_]+)"), "");
        content2.push_back(str2);
    }
    return content2;
}

std::map<token, std::vector<token>> make_word_could(const std::set<token> &set){
    std::map<token, std::vector<token>> r;
    for(auto &i : set){
        token t = i;
        t.clear_unique_data();
        r[t].push_back(i);
    }
    return std::move(r);
}

std::string make_sentence(const word_graph &w, const sentences &s, const std::string &str){
    token_vector seed = token::make_vector(mecab_tagger->parseToNode(str.c_str()));
    static std::random_device rnd_device;
    static std::mt19937 rnd(rnd_device());
    std::set<token> set;
    for(auto &i : seed){
        w.get_relational_words(set, i, 1);
    }
    std::map<token, std::vector<token>> cloud = make_word_could(set);
    std::string r;
    const sentences::frame &f = s.sentence_vec[rnd() % s.sentence_vec.size()];
    for(auto &i : f.vec){
        if(i.knockout){
            token u = i.t;
            u.clear_unique_data();
            auto iter = cloud.find(u);
            if(iter == cloud.end()){
                r += i.t.data[i.t.size - 1];
            }else{
                const auto &v = iter->second;
                r += v[rnd() % v.size()].data[i.t.size - 1];
            }
        }else{
            r += i.t.data[i.t.size - 1];
        }
    }
    return r;
}

void init_rc(){
    user_name = "asurino_bot";
    std::ofstream ofile("./asurino_rc");
    ofile << user_name;
}

void load_rc(){
    std::ifstream ifile("./asurino_rc");
    if(ifile.fail()){
        init_rc();
    }else{
        ifile >> user_name;
    }
}

int main (){
    load_rc();

    word_graph w;
    sentences s;

    {
        std::ifstream ifile("word_graph");
        if(!ifile.fail()){
            w.restore(ifile);
        }
    }
    {
        std::ifstream ifile("sentences");
        if(!ifile.fail()){
            s.restore(ifile);
        }
    }

    std::chrono::time_point<std::chrono::system_clock> prev;
    prev = std::chrono::system_clock::now() - std::chrono::seconds(1800 + 10);
    while(true){
        if(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - prev).count() >= 1800){
            prev = std::chrono::system_clock::now();

            s.sentence_vec.clear();
            std::system(("ctw -u " + user_name + " -r > log.txt").c_str());
            std::ifstream ifile("log.txt");
            auto log = read_log(ifile);
            for(auto &i : log){
                auto v = token::make_vector(mecab_tagger->parseToNode(i.c_str()));
                w.set_str(v);
                s.add_sentence(v);
            }
            std::string r = make_sentence(w, s, log.back());
            std::system(("ctw -u " + user_name + " -p \"" + r + "\"").c_str());

            {
                std::ofstream ofile("word_graph");
                w.serialize(ofile);
            }
            {
                std::ofstream ofile("sentences");
                s.serialize(ofile);
            }
        } 
        std::this_thread::sleep_for(std::chrono::seconds(1));
    }

    return 0;
}
