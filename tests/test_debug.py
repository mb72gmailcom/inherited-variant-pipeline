from inherited.debug import get_memory_usage_mb, log_memory_if_due


def test_get_memory_usage_mb_returns_positive():
    assert get_memory_usage_mb() > 0


def test_log_memory_if_due_only_when_debug_enabled(capsys):
    log_memory_if_due(50_000, debug=False, memory_block=50_000)
    assert capsys.readouterr().out == ""

    log_memory_if_due(50_000, debug=True, memory_block=50_000)
    output = capsys.readouterr().out
    assert "[debug] variants=50,000" in output
    assert "memory=" in output


def test_log_memory_if_due_respects_memory_block(capsys):
    log_memory_if_due(49_999, debug=True, memory_block=50_000)
    assert capsys.readouterr().out == ""

    log_memory_if_due(100_000, debug=True, memory_block=50_000)
    output = capsys.readouterr().out
    assert "variants=100,000" in output
