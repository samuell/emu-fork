import pytest
import subprocess
import os
import pandas as pd


@pytest.mark.parametrize(
    "emu_command, abundance_file, expected_abundancies",
    [
        (
            "./emu abundance --db emu_database --output-dir test_results --keep-files example/full_length.fa",
            "test_results/full_length_rel-abundance.tsv",
            {
                "Sphingobacterium puteale": 0.5,
                "Mycobacterium saskatchewanense": 0.25,
                "Streptococcus sobrinus": 0.25,
            },
        ),
        (
            "./emu abundance --db emu_database --output-dir test_results --type sr --keep-files example/short_read_f.fq",
            "test_results/short_read_f_rel-abundance.tsv",
            {
                "Staphylococcus hominis": 1.0,
            },
        ),
        (
            "./emu abundance --db emu_database --output-dir test_results --type sr --keep-files example/short_read_f.fq example/short_read_r.fq",
            "test_results/short_read_f-short_read_r_rel-abundance.tsv",
            {
                "Staphylococcus aureus": 1.0 / 3,
                "Salmonella enterica": 1.0 / 3,
                "Enterococcus columbae": 1.0 / 3,
            },
        ),
    ],
)
def test_cli_on_example_datasets(emu_command, abundance_file, expected_abundancies):
    if os.path.exists(abundance_file):
        os.remove(abundance_file)

    subprocess.check_output(emu_command, shell=True)

    assert os.path.exists(abundance_file)

    df = abundance_df = pd.read_csv(abundance_file, sep="\t")

    for species, expected_abundance in expected_abundancies.items():
        assert df[df["species"] == species].iloc[0]["abundance"] == expected_abundance
