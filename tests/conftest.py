import os
import sys

# Allow tests to import the package from src/ without an install step.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
