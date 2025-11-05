import sys
import subprocess

# List of required packages
required_packages = ["pandas", "gseapy", "matplotlib", "matplotlib-venn"]

for pkg in required_packages:
    try:
        __import__(pkg.replace("-", "_"))
        print(f" Package '{pkg}' is already installed.")
    except ImportError:
        print(f" Installing '{pkg}'...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])
