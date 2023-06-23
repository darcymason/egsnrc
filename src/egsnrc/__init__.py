"""EGSnrc radiation transport simulation package in Python with GPU acceleration"""
import logging.config
from pathlib import Path
from egsnrc._version import __version__, __version_info__


HERE = Path(__file__).resolve().parent

logging.config.fileConfig(HERE / "logging.conf")
