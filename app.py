import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from io import StringIO

st.set_page_config(page_title="Nucleo｜DNA 序列分析", page_icon="🧬", layout="wide")

st.markdown("""
<style>
.stApp { background: #0b1018; color: #f3f4f6; }
h1, h2, h3, p, label { color: #f3f4f6; }
div[data-testid="stMetric"] { background: #151c29; padding: 14px; border-radius: 12px; }
</style>
""", unsafe_allow_html=True)

CODON_TABLE = {
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}
AA3 = {'A':'Ala','R':'Arg','N':'Asn','D':'Asp','C':'Cys','Q':'Gln','E':'Glu','G':'Gly','H':'His','I':'Ile','L':'Leu','K':'Lys','M':'Met','F':'Phe','P':'Pro','S':'Ser','T':'Thr','W':'Trp','Y':'Tyr','V':'Val','*':'Stop'}

def clean_sequence(text):
    seq = ''.join(text.upper().split())
    invalid = sorted(set(seq) - set("ATCG"))
    if invalid:
        raise ValueError("序列中包含不支持的字符：" + "、".join(invalid) + "。目前只接受 A、T、C、G。")
    return seq

def reverse_complement(seq):
    return seq.translate(str.maketrans("ATCG", "TAGC"))[::-1]

def translate(seq, frame=0, stop_at_stop=False):
    seq = seq[frame:]
    aas, codons = [], []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, '?')
        codons.append(codon)
        aas.append(aa)
        if aa == '*' and stop_at_stop:
            break
    return codons, aas

st.title("🧬 NUCLEO")
st.subheader("DNA 序列分析工作台")
st.write("输入 DNA 序列，即可查看互补链、RNA 转录结果、氨基酸翻译和碱基组成。")

with st.sidebar:
    st.header("分析设置")
    frame_label = st.selectbox("翻译阅读框", ["阅读框 1（从第 1 个碱基开始）", "阅读框 2（从第 2 个碱基开始）", "阅读框 3（从第 3 个碱基开始）"])
    frame = ["阅读框 1（从第 1 个碱基开始）", "阅读框 2（从第 2 个碱基开始）", "阅读框 3（从第 3 个碱基开始）"].index(frame_label)
    stop_at_stop = st.toggle("遇到终止密码子时停止翻译", value=False)
    show_three = st.toggle("显示氨基酸三字母缩写", value=True)
    st.divider()
    st.header("快速示例")
    if st.button("载入示例 DNA 序列", use_container_width=True):
        st.session_state["dna_input"] = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

default_seq = st.session_state.get("dna_input", "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
raw_seq = st.text_area("DNA 序列输入框", value=default_seq, height=150, placeholder="例如：ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
st.caption("支持的字符：A、T、C、G。空格和换行会自动忽略。")

if st.button("开始分析", type="primary", use_container_width=True):
    try:
        seq = clean_sequence(raw_seq)
        if not seq:
            st.warning("请先输入一段 DNA 序列。")
            st.stop()

        comp = seq.translate(str.maketrans("ATCG", "TAGC"))
        revcomp = reverse_complement(seq)
        rna = seq.replace("T", "U")
        codons, aas = translate(seq, frame, stop_at_stop)
        protein1 = ''.join(aas)
        protein3 = ' – '.join(AA3.get(a, a) for a in aas) if show_three else None
        counts = Counter(seq)
        gc = (counts['G'] + counts['C']) / len(seq) * 100

        st.success(f"分析完成！序列长度：{len(seq)} 个碱基")
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("序列长度", f"{len(seq)} bp")
        c2.metric("GC 含量", f"{gc:.2f}%")
        c3.metric("密码子数量", len(codons))
        c4.metric("氨基酸数量", len(aas))

        st.header("1. DNA 序列结果")
        st.text_area("输入的 DNA 序列（5′ → 3′）", seq, height=80)
        st.text_area("DNA 互补链（与原链配对）", comp, height=80)
        st.text_area("DNA 反向互补链（5′ → 3′）", revcomp, height=80)

        st.header("2. RNA 转录")
        st.caption("此处按编码链规则转录：将 DNA 中的 T 替换为 U。")
        st.text_area("RNA 序列（5′ → 3′）", rna, height=80)

        st.header("3. 蛋白质翻译")
        st.caption(f"当前使用：{frame_label}。标准遗传密码表；“*”代表终止密码子。")
        st.text_area("密码子序列", ' '.join(codons), height=80)
        st.text_area("氨基酸单字母序列", protein1 or "（没有完整密码子）", height=80)
        if show_three:
            st.text_area("氨基酸三字母缩写", protein3 or "（没有完整密码子）", height=100)

        st.header("4. 碱基组成")
        bases = ['A', 'T', 'C', 'G']
        comp_df = pd.DataFrame({"碱基": bases, "数量": [counts[b] for b in bases], "比例（%）": [counts[b] / len(seq) * 100 for b in bases]})
        st.dataframe(comp_df, use_container_width=True, hide_index=True)
        fig, ax = plt.subplots()
        ax.bar(bases, [counts[b] for b in bases])
        ax.set_title("DNA 碱基数量")
        ax.set_xlabel("碱基")
        ax.set_ylabel("数量")
        fig.patch.set_alpha(0)
        ax.set_facecolor("#0b1018")
        ax.tick_params(colors="white")
        ax.xaxis.label.set_color("white")
        ax.yaxis.label.set_color("white")
        ax.title.set_color("white")
        st.pyplot(fig)
        plt.close(fig)

        report = f"""NUCLEO DNA 序列分析报告
================================
序列长度：{len(seq)} bp
GC 含量：{gc:.2f}%

输入 DNA（5'->3'）：
{seq}

DNA 互补链：
{comp}

DNA 反向互补链（5'->3'）：
{revcomp}

RNA 转录（编码链规则）：
{rna}

翻译阅读框：{frame + 1}
密码子：
{' '.join(codons)}

氨基酸单字母序列：
{protein1}

碱基计数：
A={counts['A']}, T={counts['T']}, C={counts['C']}, G={counts['G']}
"""
        st.download_button("下载分析报告（TXT）", report, file_name="nucleo_DNA分析报告.txt", mime="text/plain", use_container_width=True)
    except ValueError as e:
        st.error(str(e))

st.divider()
st.caption("Nucleo｜入门级 DNA 序列分析工具。结果用于学习参考；RNA 转录默认按编码链规则处理。")
