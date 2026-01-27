REQUIRED_PACKAGES = [
    "pandas",
    "numpy",
    "matplotlib",
    "seaborn",
    "scipy",
    "pyranges",
    "pybedtools",
    "tqdm",
]

failed = []

for pkg in REQUIRED_PACKAGES:
    try:
        __import__(pkg)
        print(f"[OK] {pkg}")
    except ImportError as e:
        print(f"[FAIL] {pkg}: {e}")
        failed.append(pkg)

if failed:
    raise SystemExit(
        f"\nMissing packages: {', '.join(failed)}\n"
        "Please install them before running the pipeline."
    )
else:
    print("\nAll required packages are installed ✅")
