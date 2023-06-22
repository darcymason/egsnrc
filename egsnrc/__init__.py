import logging.config
from pathlib import Path

HERE  = Path(__file__).resolve().parent

logging.config.fileConfig(HERE / 'logging.conf')
