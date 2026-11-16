import pathlib
import subprocess

BENCH = pathlib.Path(__file__).parent.parent / "benchmarks"


def test_no_hard_coded_paths():
    for path in BENCH.glob("*"):
        if path.suffix in (".sh", ".py"):
            text = path.read_text()
            assert "/Users/" not in text, path
            assert path.name == "config.sh" or "/Volumes/" not in text, path


def test_shell_scripts_parse():
    for path in BENCH.glob("*.sh"):
        subprocess.run(["zsh", "-n", str(path)], check=True)


def test_python_scripts_compile():
    for path in BENCH.glob("*.py"):
        compile(path.read_text(), str(path), "exec")


def test_scripts_present():
    names = {p.name for p in BENCH.glob("*")}
    for expected in ("config.sh", "bench.sh", "run_org.sh", "run_sani.sh", "split.py"):
        assert expected in names
