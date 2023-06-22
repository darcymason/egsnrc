"""EGSnrc radiation transport simulation package in Python with GPU acceleration"""
import logging.config
from pathlib import Path
import re


__version__ = "0.1.0"

__version_info__ = tuple(
    re.match(r"(\d+\.\d+\.\d+).*", __version__).group(1).split(".")
)


HERE = Path(__file__).resolve().parent

logging.config.fileConfig(HERE / "logging.conf")
