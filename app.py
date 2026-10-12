import streamlit as st
import pandas as pd
import plotly.express as px

st.set_page_config(
    page_title="Nucleo | DNA Lab",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ---------- Styling ----------
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Mono:wght@400;500&family=Manrope:wght@400;500;600;700;800&display=swap');
:root { color-scheme: dark; }
html, body, [class*="css"] { font-family: 'Manrope', sans-serif; }
.stApp {
    background:
      radial-gradient(ellipse at 12% 0%, rgba(38, 112, 100, .18), transparent 36%),
      radial-gradient(ellipse at 92% 8%, rgba(91, 76, 177, .16), transparent 30%),
      #090d13;
    color: #e8edf5;
}
.block-container { padding-top: 2.2rem; max-width: 1440px; }
[data-testid="stSidebar"] { background: #0c121b; border-right: 1px solid #1c2735; }
.hero-kicker { color: #6ee7c1; font: 500 12px 'DM Mono', monospace; letter-spacing: .2em; text-transform: uppercase; }
.hero-title { font-size: clamp(2.4rem, 5vw, 4.2rem); font-weight: 800; letter-spacing: -.055em; line-height: 1.02; margin: .45rem 0 .8rem; color: #f4f7fb; }
.hero-title span { color: #6ee7c1; }
.hero-sub { color: #96a4b8; font-size: 1rem; max-width: 680px; line-height: 1.7; }
.section-label { color: #a9b8cb; font: 500 11px 'DM Mono', monospace; letter-spacing: .16em; text-transform: uppercase; margin: .5rem 0 .65rem; }
.metric-card { background: linear-gradient(145deg, #131d29, #0d141e); border: 1px solid #243344; border-radius: 16px; padding: 18px 20px; min-height: 108px; }
.metric-label { color: #8ea0b5; font-size: 12px; text-transform: uppercase; letter-spacing: .1em; }
.metric-value { color: #f1f6fb; font-size: 2rem; font-weight: 800; letter-spacing: -.04em; margin-top: 5px; }
.metric-note { color: #6ee7c1; font-size: 11px; margin-top: 2px; }
.panel { background: #0e151f; border: 1px solid #243344; border-radius: 16px; padding: 18px 20px; }
.panel-title { color: #edf3fa; font-size: 1.03rem; font-weight: 700; margin-bottom: .3rem; }
.panel-caption { color: #8495aa; font-size: .82rem; margin-bottom: 1rem; }
.sequence-box { background: #080d13; border: 1px solid #1f3340; border-radius: 12px; padding: 14px; color: #80f0c9; font: 13px/1.9 'DM Mono', monospace; overflow-wrap: anywhere; }
.small-muted { color: #7f90a5; font-size: .82rem; }
div.stButton > button { border-radius: 10px; font-weight: 700; min-height: 2.8rem; }
.stTextArea textarea { background: #080d13 !important; color: #dce7f4 !important; border-color: #263748 !important; font-family: 'DM Mono', monospace !important; }
.stTextInput input { background: #080d13 !important; color: #dce7f4 !important; border-color: #263748 !important; }
hr { border-color: #202d3b; }
</style>
""", unsafe_allow_html=True)

# ---------- Biology logic ----------
DNA_COMPLEMENT = str.maketrans({"A": "T", "T": "A", "C": "G", "G": "C"})

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
AA3 = {
    'A':'Ala','R':'Arg','N':'Asn','D':'Asp','C':'Cys','Q':'Gln','E':'Glu','G':'Gly',
    'H':'His','I':'Ile','L':'Leu','K':'Lys','M':'Met','F':'Phe','P':'Pro','S':'Ser',
    'T':'Thr','W':'Trp','Y':'Tyr','V':'Val','*':'Stop'
}

def translate_sequence(dna, frame=0, stop_at_stop=False):
    usable = dna[frame:]
    aa = []
    codons = []
    for i in range(0, len(usable) - 2, 3):
        codon = usable[i:i+3]
        codons.append(codon)
        residue = CODON_TABLE.get(codon, '?')
        aa.append(residue)
        if residue == '*' and stop_at_stop:
            break
    return codons, ''.join(aa)

def pretty_codons(codons):
    return ' '.join(codons) if codons else '—'

# ---------- Sidebar ----------
with st.sidebar:
    st.markdown("## 🧬 NUCLEO")
    st.caption("SEQUENCE ANALYSIS WORKBENCH")
    st.markdown("---")
    st.markdown("### Analysis settings")
    frame_label = st.selectbox("Translation reading frame", ["Frame 1 (start at base 1)", "Frame 2 (start at base 2)", "Frame 3 (start at base 3)"])
    frame = ["Frame 1 (start at base 1)", "Frame 2 (start at base 2)", "Frame 3 (start at base 3)"].index(frame_label)
    stop_at_stop = st.toggle("Stop translation at stop codon", value=False)
    show_three_letter = st.toggle("Show three-letter amino-acid codes", value=True)
    st.markdown("---")
    st.markdown("### Quick sample")
    if st.button("Load sample sequence", use_container_width=True):
        st.session_state["dna_input"] = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    st.markdown("<div class='small-muted'>Accepted alphabet: A, T, C, G. Spaces and line breaks are ignored.</div>", unsafe_allow_html=True)

# ---------- Header ----------
st.markdown("<div class='hero-kicker'>MOLECULAR BIOLOGY · TOOLKIT 001</div>", unsafe_allow_html=True)
st.markdown("<div class='hero-title'>DNA, decoded<span>.</span></div>", unsafe_allow_html=True)
st.markdown("<div class='hero-sub'>A compact sequence workbench for complementary strands, transcription, translation, and GC composition. Paste a sequence and explore its molecular readout.</div>", unsafe_allow_html=True)
st.write("")

left, right = st.columns([1.08, .92], gap="large")
with left:
    st.markdown("<div class='section-label'>01 / Input sequence</div>", unsafe_allow_html=True)
    default_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    raw = st.text_area(
        "DNA sequence",
        value=st.session_state.get("dna_input", default_seq),
        height=150,
        placeholder="Paste DNA sequence here…",
        label_visibility="collapsed",
        help="Paste a DNA sequence containing A, T, C, and G. Whitespace is ignored."
    )
    cleaned = ''.join(raw.upper().split())
    invalid = sorted(set(cleaned) - set("ATCG"))
    run = st.button("↗  Analyse sequence", type="primary", use_container_width=True)
    if run:
        st.session_state["last_sequence"] = cleaned

with right:
    st.markdown("<div class='section-label'>02 / Sequence health</div>", unsafe_allow_html=True)
    if not cleaned:
        st.info("Enter a DNA sequence to begin analysis.")
    elif invalid:
        st.error("Invalid characters detected: " + ", ".join(invalid) + ". Use only A, T, C, G.")
    else:
        gc = (cleaned.count("G") + cleaned.count("C")) / len(cleaned) * 100
        base_df = pd.DataFrame({"Base": ["A", "T", "G", "C"], "Count": [cleaned.count(b) for b in "ATGC"]})
        c1, c2 = st.columns(2)
        with c1:
            st.markdown(f"<div class='metric-card'><div class='metric-label'>Sequence length</div><div class='metric-value'>{len(cleaned):,}</div><div class='metric-note'>nucleotides</div></div>", unsafe_allow_html=True)
        with c2:
            st.markdown(f"<div class='metric-card'><div class='metric-label'>GC content</div><div class='metric-value'>{gc:.1f}%</div><div class='metric-note'>G + C / total bases</div></div>", unsafe_allow_html=True)
        st.write("")
        st.progress(min(max(gc / 100, 0.0), 1.0), text=f"GC composition · {gc:.2f}%")
        st.plotly_chart(
            px.bar(base_df, x="Base", y="Count", text="Count", template="plotly_dark",
                   color="Base", color_discrete_map={"A":"#f5b84b","T":"#e77b91","G":"#6ee7c1","C":"#8f9cff"})
              .update_layout(height=210, margin=dict(l=8,r=8,t=10,b=8), paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)", showlegend=False,
                             xaxis_title=None, yaxis_title=None, font=dict(family="DM Mono", color="#aebdce"))
              .update_traces(textposition="outside", cliponaxis=False),
            use_container_width=True
        )

st.markdown("---")

if not cleaned:
    st.stop()
if invalid:
    st.warning("Fix the invalid characters above to see the sequence outputs.")
    st.stop()

comp = cleaned.translate(DNA_COMPLEMENT)
rev_comp = comp[::-1]
rna = cleaned.replace("T", "U")
codons, aa = translate_sequence(cleaned, frame=frame, stop_at_stop=stop_at_stop)
aa3 = ' – '.join(AA3.get(x, '?') for x in aa) if aa else '—'
complete_codons = len(cleaned[frame:]) // 3
remainder = len(cleaned[frame:]) % 3

st.markdown("<div class='section-label'>03 / Molecular outputs</div>", unsafe_allow_html=True)
tab1, tab2, tab3, tab4 = st.tabs(["Complementary strand", "Reverse complement", "RNA transcript", "Protein translation"])

def sequence_panel(title, subtitle, sequence, badge):
    st.markdown(f"<div class='panel'><div class='panel-title'>{title}</div><div class='panel-caption'>{subtitle}</div><div class='sequence-box'>{sequence}</div><div class='small-muted' style='margin-top:10px'>{badge}</div></div>", unsafe_allow_html=True)

with tab1:
    sequence_panel("Complementary DNA strand", "Base pairing: A ↔ T · C ↔ G", comp, f"{len(comp):,} bases")
with tab2:
    sequence_panel("Reverse-complement DNA", "Complement, then reverse the strand direction (5′ → 3′)", rev_comp, f"{len(rev_comp):,} bases")
with tab3:
    sequence_panel("RNA transcript (coding-strand convention)", "T is replaced by U; assumes the input is the coding DNA strand", rna, f"{len(rna):,} nucleotides")
with tab4:
    st.markdown("<div class='panel-title'>Translation</div>", unsafe_allow_html=True)
    st.markdown(f"<div class='panel-caption'>Reading frame: {frame + 1} · {complete_codons} complete codons · {remainder} trailing base(s) excluded</div>", unsafe_allow_html=True)
    st.markdown("<div class='section-label'>Codons</div>", unsafe_allow_html=True)
    st.markdown(f"<div class='sequence-box'>{pretty_codons(codons)}</div>", unsafe_allow_html=True)
    st.markdown("<div class='section-label' style='margin-top:16px'>Amino-acid sequence · one-letter</div>", unsafe_allow_html=True)
    st.markdown(f"<div class='sequence-box'>{aa or '—'}</div>", unsafe_allow_html=True)
    if show_three_letter:
        st.markdown("<div class='section-label' style='margin-top:16px'>Amino-acid sequence · three-letter</div>", unsafe_allow_html=True)
        st.markdown(f"<div class='sequence-box'>{aa3}</div>", unsafe_allow_html=True)
    st.caption("Translation uses the standard genetic code. * marks a stop codon; ? marks an unknown codon.")

st.markdown("---")
st.markdown("<div class='section-label'>04 / Base composition</div>", unsafe_allow_html=True)
base_counts = {b: cleaned.count(b) for b in "ATGC"}
total = len(cleaned)
cols = st.columns(4)
for col, base in zip(cols, "ATGC"):
    pct = base_counts[base] / total * 100
    with col:
        st.markdown(f"<div class='metric-card'><div class='metric-label'>Base {base}</div><div class='metric-value'>{base_counts[base]:,}</div><div class='metric-note'>{pct:.1f}% of sequence</div></div>", unsafe_allow_html=True)

st.write("")
st.download_button(
    "↓  Download analysis as TXT",
    data=(
        f"NUCLEO DNA ANALYSIS\nLength: {len(cleaned)} bp\nGC content: {(cleaned.count('G')+cleaned.count('C'))/len(cleaned)*100:.2f}%\n\n"
        f"Input DNA:\n{cleaned}\n\nComplement:\n{comp}\n\nReverse complement:\n{rev_comp}\n\nRNA (coding-strand convention):\n{rna}\n\n"
        f"Codons (frame {frame+1}):\n{pretty_codons(codons)}\n\nProtein (one-letter):\n{aa}\n\nProtein (three-letter):\n{aa3}\n"
    ),
    file_name="nucleo_analysis.txt",
    mime="text/plain",
)

st.caption("NUCLEO / Educational use · Translation assumes the standard genetic code and the selected reading frame. RNA output assumes the input DNA is the coding strand.")
