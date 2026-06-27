from inherited.output import BlockWriter, serialize_patient_ids


def test_block_writer_flushes_in_blocks(tmp_path):
    path = tmp_path / "out.tsv"
    writer = BlockWriter(path, block_size=2)
    hits = {"p1": ("0/1", "0/0", "0/1", "30")}

    writer.append("22", "100", "A", "G", serialize_patient_ids(hits))
    writer.append("22", "200", "C", "T", serialize_patient_ids(hits))
    assert writer.lines_written == 2
    assert writer.block == []

    writer.append("22", "300", "G", "A", serialize_patient_ids(hits))
    assert writer.lines_written == 2
    assert len(writer.block) == 1

    writer.close()
    lines = path.read_text(encoding="utf-8").strip().splitlines()
    assert lines[0].startswith("#CHROM")
    assert len(lines) == 4
