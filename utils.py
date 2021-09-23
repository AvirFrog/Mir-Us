"""
Utility functions used in other modules.

"""

import functools
import os
import datetime
import traceback
from timeit import default_timer as timer

import dill
from colorama import init, Fore

__authors__ = ["Kacper Dudczak, Maciej Michalczyk"]
__copyright__ = "Copyright 2021, mirBase Project"
__credits__ = ["Kacper Dudczak", "Maciej Michalczyk", "Marta Wysocka", "Marek Żywicki"]
__license__ = "MIT"
__version__ = "0.5"
__maintainer__ = ["Kacper Dudczak", "Maciej Michalczyk"]
__email__ = ["kacper.dudczak19@gmail.com", "mccv99@gmail.com"]
__status__ = "Production"
__deprecated__ = False

# GLOBALS---------------------------------------------------------------------
std_banner = f'''

        {Fore.LIGHTGREEN_EX}0           o 0           o 0           o {Fore.LIGHTGREEN_EX}ooo        ooooo  o8o                   ooooo     ooo
        {Fore.GREEN}| 0       o | | 0       o | | 0       o | {Fore.GREEN}`88.       .888'  `"'                   `888'     `8'
        {Fore.YELLOW}| | 0   o | | | | 0   o | | | | 0   o | | {Fore.YELLOW} 888b     d'888  oooo  oooo d8b          888       8   .oooo.o
        {Fore.RED}| | | 0 | | | | | | 0 | | | | | | 0 | | | {Fore.RED} 8 Y88. .P  888  `888  `888""8P          888       8  d88(  "8
        {Fore.LIGHTRED_EX}| | o   0 | | | | o   0 | | | | o   0 | | {Fore.LIGHTRED_EX} 8  `888'   888   888   888     8888888  888       8  `"Y88b.
        {Fore.MAGENTA}| o       0 | | o       0 | | o       0 | {Fore.MAGENTA} 8    Y     888   888   888              `88.    .8'  o.  )88b
        {Fore.CYAN}o           0 o           0 o           0 {Fore.CYAN}o8o        o888o o888o d888b               `YbodP'    8""888P'
        '''


# UTILITY FUNCTIONS-----------------------------------------------------------
def time_this(func):
    """
    Decorator which returns information about execution of decorated function.

    :param func: Any miBase function
    :return: Execution time and values returned by a function
    """

    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start = timer()
        values = func(*args, **kwargs)
        end = timer()
        runtime = end - start
        if values is None:
            print(f"{Fore.RED}[Mir-Us]  {func.__name__!r} No records matching given criteria.")
        else:
            print(f"{Fore.GREEN}[Mir-Us]  {func.__name__!r} found {values[1]} results in {runtime:.6f} seconds")
            return values[0]

    return wrapper_timer


def _fatal_error_handle(e):
    msg = f"\n{Fore.RED}[Mir-Us]   Error! Cannot compile data files:"
    if isinstance(e, KeyError):
        msg = msg + f"{Fore.RED} Given database version is not compatible with the tool or unexpected " \
                    f"key appeared.\n"
    os.makedirs(f"{os.getcwd()}/logs", exist_ok=True)
    time = datetime.datetime.now().strftime("%m_%d_%Y-%H_%M_%S")
    with open(f"logs/error_{time}.log", "a") as f_log:
        f_log.write(traceback.format_exc())
    return msg


def _cache_versions():
    """Caches metadata about miRBase database versions

    Returns:
        dict: miRBase database metadata
    """
    with open("versions.mir", 'rb') as ver:
        return dill.loads(ver.read())


def _exists(to_compare, obj):
    """Utility function, which compares list of objects to particular object, to check if compared object is present
     in the list

    Args:
        to_compare (list[Any]): List of objects to compare
        obj (Any): Object to be compared with list of objects

    Returns:
        bool: True if compared object is present in the list of objects
    """
    for elem in to_compare:
        if elem is obj:
            return True


def _show_banner():
    """
    Utility function, which displays Mir-Us banner at start.
    """
    print(std_banner)
