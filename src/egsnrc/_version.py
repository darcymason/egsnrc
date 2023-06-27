import re
from typing import cast, Match
from importlib.metadata import version

__version__: str = version("egsnrc")

result = cast(Match[str], re.match(r'(\d+\.\d+\.\d+).*', __version__))
__version_info__ = tuple(result.group(1).split('.'))
