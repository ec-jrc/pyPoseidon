import pathlib

# specify package paths
ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent
TEST_DIR = ROOT_DIR / "tests"
DATA_DIR = TEST_DIR / "data"
ROTATED_DIR = DATA_DIR / "rotated"

DATA_DIR.mkdir(exist_ok=True)
ROTATED_DIR.mkdir(exist_ok=True)
