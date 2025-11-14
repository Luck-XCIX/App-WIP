import streamlit as st
import tempfile
import os
from primer_designer import create_design
from datetime import datetime

def run_primer_design(fasta_path, primer_type, design_mode, qpcr_params, race_params):
    return create_design(fasta_path, primer_type, design_mode.lower(), qpcr_params, race_params)


# -------------------------------
# Session State Initialization
# -------------------------------
if "output_df1" not in st.session_state:
    st.session_state.output_df1 = None
if "output_df2" not in st.session_state:
    st.session_state.output_df2 = None
if "fasta_path" not in st.session_state:
    st.session_state.fasta_path = None


# -------------------------------
# Clear output when input changes
# -------------------------------
def clear_output():
    st.session_state.output_df1 = None
    st.session_state.output_df2 = None


# -------------------------------
# UI
# -------------------------------
st.title("🧬 PRIMEGEN")

st.markdown("""
PRIMEGEN is a Python-based primer design tool that identifies optimal primer pairs from multi-isoform FASTA sequences for **qPCR** or **RACE** applications. It operates in two distinct modes:

* **Specificity Mode**: Identifies primer pairs that are unique to a single target isoform, ensuring minimal cross-reactivity. The final output for this mode includes a detailed analysis of primer quality, including GC content, melting temperature (Tm), and predicted secondary structures (hairpins, self-dimers, and cross-dimers) to assess biological viability.

* **Coverage Mode**: Addresses the challenge of amplifying all isoforms with a minimum number of primers. It employs a greedy algorithm to solve this "Set Cover" type problem, prioritizing primer pairs that cover the largest number of previously unamplified transcripts in each step.

The pipeline filters all potential primers based on user-defined parameters for GC content, melting temperature, and k-mer length. For qPCR, it also enforces a specific product size range. Its use of cached functions for primer hit detection ensures efficient and scalable analysis, even with large datasets.
""")

st.divider()

uploaded_file = st.file_uploader(
    "📁 Choose FASTA file", type=["fa", "fasta"], on_change=clear_output
)
primer_type = st.selectbox(
    "🧪 Select Primer Type", ["RACE", "qPCR"], on_change=clear_output
)

race_params = {}
if primer_type == "RACE":
    with st.expander("⚙️ Change design parameters", expanded=False):
        col1, col2 = st.columns(2, gap="large")
        race_params["window_size"] = col1.number_input("🔍 Search Window Size (bp)", value=200, min_value=50, on_change=clear_output)
        race_params["k_range"] = col1.slider("📏 Primer Length (k-mer)", value=(23, 29), min_value=15, max_value=30, on_change=clear_output)
        race_params["gc_range"] = col2.slider("🧬 GC (%)", value=(0.5,0.7), min_value=0.0, max_value=1.0, on_change=clear_output)
        race_params["tm_range"] = col2.slider("🌡️ Tm (°C)", value=(57.0,62.0), step=0.1, min_value=50.0, max_value=65.0, on_change=clear_output)

qpcr_params = {}
if primer_type == "qPCR":
    with st.expander("⚙️ Change design parameters", expanded=False):
        col1, col2 = st.columns(2, gap="large")
        qpcr_params["fw_range"] = col1.slider("⏩ Forward Search Region", value=(50,100), min_value=50, max_value=150, on_change=clear_output)
        qpcr_params["rev_range"] = col1.slider("⏪ Reverse Search Region", value=(150,200), min_value=150, max_value=250, on_change=clear_output)
        qpcr_params["k_range"] = col1.slider("📏 Primer Length (k-mer)", value=(18,26), min_value=15, max_value=30, on_change=clear_output)
        qpcr_params["gc_range"] = col2.slider("🧬 GC (%)", value=(0.5,0.7), min_value=0.0, max_value=1.0, on_change=clear_output)
        qpcr_params["prod_range"] = col2.slider("💥 Product Size (bp)", value=(80,150), min_value=50, max_value=200, on_change=clear_output)
        qpcr_params["tm_range"] = col2.slider("🌡️ Tm (°C)", value=(57.0,62.0), step=0.1, min_value=50.0, max_value=65.0, on_change=clear_output)

design_mode = st.selectbox(
    "🎨 Select Design Mode", ["Specificity", "Coverage"], on_change=clear_output
)

# -------------------------------
# Run Primer Design
# -------------------------------
if st.button("⚡ Run Primer Design", use_container_width=True):
    if uploaded_file is None:
        st.error("Please choose a FASTA file first.")
    else:
        start_time = datetime.now()
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fa") as tmp_file:
            tmp_file.write(uploaded_file.read())
            st.session_state.fasta_path = tmp_file.name

        with st.spinner("Running primer design..."):
            num_sequences, output_df1, output_df2 = run_primer_design(
                st.session_state.fasta_path, 
                primer_type, 
                design_mode, 
                qpcr_params=qpcr_params if primer_type == "qPCR" else None,
                race_params=race_params if primer_type == "RACE" else None,
            )

        duration = datetime.now() - start_time
        duration = duration.total_seconds()
        seq_per_s = int(round(num_sequences / duration, 0))
        duration *= 1000

        if duration < 1000:
            duration_text = f"{int(duration)} ms"
        elif duration < 60000:
            duration_text = f"{duration / 1000:.2f} seconds"
        else:
            duration_text = f"{duration / 60000:.2f} minutes"

        st.divider()

        if output_df1 is None:
            st.error(f"No suitable primers were found after {duration_text}.")
        else:
            st.success(f"Primer design completed in {duration_text} ({num_sequences} sequences at {seq_per_s} seq/s)")
            st.session_state.output_df1 = output_df1
            st.session_state.output_df2 = output_df2

# -------------------------------
# Display Results
# -------------------------------

if st.session_state.output_df1 is not None:
    pre = "Best" if design_mode == "Specificity" else "Selected"
    num_pairs = st.session_state.output_df1.shape[0]
    pairs = f"pair{'s' if num_pairs != 1 else ''}"
    st.subheader(f"{pre} design ({num_pairs} {pairs})")
    st.dataframe(st.session_state.output_df1)

    csv1 = st.session_state.output_df1.to_csv(index=False).encode("utf-8")
    st.download_button(
        label=f"⬇️ Download {pre.lower()} design",
        data=csv1,
        file_name=f"{pre}_design_{primer_type}_{design_mode}.csv".lower(),
        mime="text/csv",
    )

if st.session_state.output_df2 is not None:
    num_pairs = st.session_state.output_df2.shape[0]
    pairs = f"pair{'s' if num_pairs != 1 else ''}"
    st.subheader(f"Others ({num_pairs} {pairs})")
    st.dataframe(st.session_state.output_df2)

    csv2 = st.session_state.output_df2.to_csv(index=False).encode("utf-8")
    st.download_button(
        label=f"⬇️ Download others",
        data=csv2,
        file_name=f"others_{primer_type}_{design_mode}.csv".lower(),
        mime="text/csv",
    )

# -------------------------------
# Clean Up Temp File
# -------------------------------
if st.session_state.fasta_path and os.path.exists(st.session_state.fasta_path):
    os.remove(st.session_state.fasta_path)
    st.session_state.fasta_path = None
