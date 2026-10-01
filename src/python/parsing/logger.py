import logging
import sys
import os

LOGGING_MAP = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
}

# logging definition
logger = logging.getLogger(__name__)
# logger.setLevel(logging.DEBUG)
formater = "%(asctime)s - %(module)s - %(levelname)s - %(message)s"
# añadimos un handler para stdout
c_handler = logging.StreamHandler(sys.stdout)
c_format = logging.Formatter(formater)
c_handler.setFormatter(c_format)
logger.addHandler(c_handler)


def file_logging(
    path: str = "_log/", filename: str = "generic", save_in_log: bool = True
):
    """add file handler
    Save in log means the file will be saved in the _log folder directly
    (_log/basename(filename).log), if False, the file will be saved in
    the same path as the output file"""
    if save_in_log:
        basen = os.path.basename(filename)
        final_path = os.path.join(path, f"{os.path.splitext(basen)[0]}.log")
    else:
        final_path = os.path.join(path, filename)
    if not os.path.exists(final_path):
        print(path)
        os.makedirs(path, exist_ok=True)
    f_handler = logging.FileHandler(final_path, mode="w")
    f_handler.setFormatter(c_format)
    logger.addHandler(f_handler)
    logger.info(" ".join(sys.argv))
