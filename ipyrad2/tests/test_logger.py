

from loguru import logger
from ipyrad2.utils.logger import set_log_level

set_log_level("DEBUG")
logger.info("THIS IS A TEST.")

set_log_level("DEBUG")
logger.debug("HI")
