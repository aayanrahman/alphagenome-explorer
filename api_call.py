import os
from dotenv import load_dotenv
from alphagenome.models.dna_client import create, OutputType
from alphagenome.data.genome import Interval, Variant

load_dotenv()

# BIN1 rs6733839: chr2:127609280, C > T
# Interval centered on the variant (500kb window)
BIN1_VARIANT = Variant(
    chromosome="chr2",
    position=127609280,
    reference_bases="C",
    alternate_bases="T",
)
BIN1_INTERVAL = Interval(
    chromosome="chr2",
    start=127609280 - 262144,
    end=127609280 + 262144,
)


def run_bin1_variant():
    api_key = os.environ["ALPHAGENOME_API_KEY"]
    client = create(api_key=api_key)

    variant_output = client.predict_variant(
        interval=BIN1_INTERVAL,
        variant=BIN1_VARIANT,
        requested_outputs=[OutputType.RNA_SEQ, OutputType.DNASE],
        # CL:0000129 (microglia) not supported by API; using frontal cortex instead
        ontology_terms=["UBERON:0000955", "UBERON:0001870"],
    )

    return variant_output


if __name__ == "__main__":
    output = run_bin1_variant()
    print(f"REF rna_seq shape:  {output.reference.rna_seq.values.shape}")
    print(f"REF dnase shape:    {output.reference.dnase.values.shape}")
    print(f"ALT rna_seq shape:  {output.alternate.rna_seq.values.shape}")
    print(f"ALT dnase shape:    {output.alternate.dnase.values.shape}")
