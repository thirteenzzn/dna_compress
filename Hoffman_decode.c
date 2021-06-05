//解码Huffman编码
void HuffmanDecoding(string s) {
    vector<string> v;
    // 标识位
    int ok = 1;
    for(int i = 0; i < s.length();) {
        // 根节点
        int x = (2 * n)-1 -1;
        // 不为叶子节点
        while(huffman[x].lchild != -1 && huffman[x].rchild != -1) {
            // 左子树
            if(s[i] == '0') {
                x = huffman[x].lchild;
            } else {
                // 右子树
                x = huffman[x].rchild;
            }
            i++;
            // 处理0,1序列有误
            // 这种情况一般是结尾0,1序列少了，导致最后一个字符串解码失败
            if(i == s.length() && huffman[x].lchild != -1) {
                ok = 0;
                break;
            }
        }

        if(ok) {
            v.push_back(huffman[x].value);
        }
    }
    if(ok) {
        for(int i = 0; i < v.size(); i++) {
            cout << v[i];
        }
        cout << endl << endl;
    } else {
        cout << "解码有误。" << endl << endl;
    }
}