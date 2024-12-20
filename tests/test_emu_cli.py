import pytest
import subprocess
import os
import pandas as pd
import tempfile


@pytest.mark.parametrize(
    "emu_command_tpl, abundance_filepath_tpl, expected_abundancies",
    [
        (
            "./emu abundance --db emu_database --output-dir {tempdir_path} --keep-files example/full_length.fa",
            "{tempdir_path}/full_length_rel-abundance.tsv",
            {
                "Sphingobacterium puteale": 0.5,
                "Mycobacterium saskatchewanense": 0.25,
                "Streptococcus sobrinus": 0.25,
            },
        ),
        (
            "./emu abundance --db emu_database --output-dir {tempdir_path} --type sr --keep-files example/short_read_f.fq",
            "{tempdir_path}/short_read_f_rel-abundance.tsv",
            {
                "Staphylococcus hominis": 1.0,
            },
        ),
        (
            "./emu abundance --db emu_database --output-dir {tempdir_path} --type sr --keep-files example/short_read_f.fq example/short_read_r.fq",
            "{tempdir_path}/short_read_f-short_read_r_rel-abundance.tsv",
            {
                "Staphylococcus aureus": 1.0 / 3,
                "Salmonella enterica": 1.0 / 3,
                "Enterococcus columbae": 1.0 / 3,
            },
        ),
    ],
)
def test_cli_on_example_datasets(emu_command_tpl, abundance_filepath_tpl, expected_abundancies):
    tempdir = tempfile.TemporaryDirectory()

    emu_command = emu_command_tpl.format(tempdir_path=tempdir.name)
    abundance_filepath = abundance_filepath_tpl.format(tempdir_path=tempdir.name)

    if os.path.exists(abundance_filepath):
        os.remove(abundance_filepath)

    subprocess.check_output(emu_command, shell=True)

    assert os.path.exists(abundance_filepath)

    df = abundance_df = pd.read_csv(abundance_filepath, sep="\t")

    for species, expected_abundance in expected_abundancies.items():
        assert df[df["species"] == species].iloc[0]["abundance"] == expected_abundance
