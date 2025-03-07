## BASIC PYTHON LIBRARIES
from subprocess import run
import textwrap
from os import path as p
##  CLEAN CODE
from typing import Optional, Union, List

from UtilitiesCloset import drListInitiator
###########################################################################################

def print_drMD_logo() -> None:
    """
    Prints the DRMD logo.

    Returns:
        None
    """
    run(["clear"])
    tealColor = "\033[38;5;37m" 
    boldText = "\033[1m"
    resetTextColor = "\033[0m"
    print(tealColor+boldText+
          """
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕
            dddddddd                                                                       
            d::::::d                   MMMMMMMM               MMMMMMMMDDDDDDDDDDDDD         
            d::::::d                   M:::::::M             M:::::::MD::::::::::::DDD     
            d::::::d                   M::::::::M           M::::::::MD:::::::::::::::DD   
            d:::::d                    M:::::::::M         M:::::::::MDDD:::::DDDDD:::::D  
    ddddddddd:::::drrrrr   rrrrrrrrr   M::::::::::M       M::::::::::M  D:::::D    D:::::D 
  dd::::::::::::::dr::::rrr:::::::::r  M:::::::::::M     M:::::::::::M  D:::::D     D:::::D
 d::::::::::::::::dr:::::::::::::::::r M:::::::M::::M   M::::M:::::::M  D:::::D     D:::::D
d:::::::ddddd:::::drr::::::rrrrr::::::rM::::::M M::::M M::::M M::::::M  D:::::D     D:::::D
d::::::d    d:::::d r:::::r     r:::::rM::::::M  M::::M::::M  M::::::M  D:::::D     D:::::D
d:::::d     d:::::d r:::::r     rrrrrrrM::::::M   M:::::::M   M::::::M  D:::::D     D:::::D
d:::::d     d:::::d r:::::r            M::::::M    M:::::M    M::::::M  D:::::D     D:::::D
d:::::d     d:::::d r:::::r            M::::::M     MMMMM     M::::::M  D:::::D    D:::::D 
d::::::ddddd::::::ddr:::::r            M::::::M               M::::::MDDD:::::DDDDD:::::D  
 d:::::::::::::::::dr:::::r            M::::::M               M::::::MD:::::::::::::::DD        
  d:::::::::ddd::::dr:::::r            M::::::M               M::::::MD::::::::::::DDD     
   ddddddddd   dddddrrrrrrr            MMMMMMMM               MMMMMMMMDDDDDDDDDDDDD

⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕
                    Molecular Dynamics: Just what the Doctor Ordered!
    """
    +resetTextColor)
###########################################################################################
def print_botched(simulationReport: List[Union[None, dict]]) -> None:
    greenText = "\033[32m"
    redText = "\033[31m"
    orangeText = "\033[38;5;172m"
    yellowText = "\033[33m"
    resetTextColor = "\033[0m"
    tealColor = "\033[38;5;37m" 

    # run(["clear"])
    print(redText+
          f"""
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕
  ██████  ██▓ ███▄ ▄███▓ █    ██  ██▓    ▄▄▄     ▄▄▄█████▓ ██▓ ▒█████   ███▄    █   ██████ 
▒██    ▒ ▓██▒▓██▒▀█▀ ██▒ ██  ▓██▒▓██▒   ▒████▄   ▓  ██▒ ▓▒▓██▒▒██▒  ██▒ ██ ▀█   █ ▒██    ▒ 
░ ▓██▄   ▒██▒▓██    ▓██░▓██  ▒██░▒██░   ▒██  ▀█▄ ▒ ▓██░ ▒░▒██▒▒██░  ██▒▓██  ▀█ ██▒░ ▓██▄   
  ▒   ██▒░██░▒██    ▒██ ▓▓█  ░██░▒██░   ░██▄▄▄▄██░ ▓██▓ ░ ░██░▒██   ██░▓██▒  ▐▌██▒  ▒   ██▒
▒██████▒▒░██░▒██▒   ░██▒▒▒█████▓ ░██████▒▓█   ▓██▒ ▒██▒ ░ ░██░░ ████▓▒░▒██░   ▓██░▒██████▒▒
▒ ▒▓▒ ▒ ░░▓  ░ ▒░   ░  ░░▒▓▒ ▒ ▒ ░ ▒░▓  ░▒▒   ▓▒█░ ▒ ░░   ░▓  ░ ▒░▒░▒░ ░ ▒░   ▒ ▒ ▒ ▒▓▒ ▒ ░
░ ░▒  ░ ░ ▒ ░░  ░      ░░░▒░ ░ ░ ░ ░ ▒  ░ ▒   ▒▒ ░   ░     ▒ ░  ░ ▒ ▒░ ░ ░░   ░ ▒░░ ░▒  ░ ░
░  ░  ░   ▒ ░░      ░    ░░░ ░ ░   ░ ░    ░   ▒    ░       ▒ ░░ ░ ░ ▒     ░   ░ ░ ░  ░  ░  
      ░   ░         ░      ░         ░  ░     ░  ░         ░      ░ ░           ░       ░  

                ▄▄▄▄    ▒█████  ▄▄▄█████▓ ▄████▄   ██░ ██ ▓█████ ▓█████▄  ▐██▌ 
                ▓█████▄ ▒██▒  ██▒▓  ██▒ ▓▒▒██▀ ▀█  ▓██░ ██▒▓█   ▀ ▒██▀ ██▌ ▐██▌ 
                ▒██▒ ▄██▒██░  ██▒▒ ▓██░ ▒░▒▓█    ▄ ▒██▀▀██░▒███   ░██   █▌ ▐██▌ 
                ▒██░█▀  ▒██   ██░░ ▓██▓ ░ ▒▓▓▄ ▄██▒░▓█ ░██ ▒▓█  ▄ ░▓█▄   ▌ ▓██▒ 
                ░▓█  ▀█▓░ ████▓▒░  ▒██▒ ░ ▒ ▓███▀ ░░▓█▒░██▓░▒████▒░▒████▓  ▒▄▄  
                ░▒▓███▀▒░ ▒░▒░▒░   ▒ ░░   ░ ░▒ ▒  ░ ▒ ░░▒░▒░░ ▒░ ░ ▒▒▓  ▒  ░▀▀▒ 
                ▒░▒   ░   ░ ▒ ▒░     ░      ░  ▒    ▒ ░▒░ ░ ░ ░  ░ ░ ▒  ▒  ░  ░ 
                ░    ░ ░ ░ ░ ▒    ░      ░         ░  ░░ ░   ░    ░ ░  ░     ░ 
                ░          ░ ░           ░ ░       ░  ░  ░   ░  ░   ░     ░    
                    ░                   ░                        ░          
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕
          """+resetTextColor)
    
    botchedSimulations = [sim for sim in simulationReport if sim["errorMessage"] is not None]
    print(f"-->{' '*4}drMD failed to complete simulations for {redText}{str(len(botchedSimulations))}{resetTextColor} out of {str(len(simulationReport))} input systems")
    print(f"-->{' '*4}Simluations on the following systems failed to complete: ")


    for botchedSimulation in botchedSimulations:
        if botchedSimulation is not None:
            print(f"{redText}{' '*7}⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕{resetTextColor}")
            print(f"{' '*7}For System:\t\t{yellowText}{botchedSimulation['pdbName']}{resetTextColor}")
            print("")
            print(f"{tealColor}{' '*7}{'#'*4}{' '*7}Traceback{' '*7}{'#'*4}{resetTextColor}")
            print(f"{' '*7}In Script:\t\t{orangeText}{botchedSimulation['scriptName']}{resetTextColor}")
            print(f"{' '*7}In Function:\t\t{orangeText}{botchedSimulation['functionName']}{resetTextColor}")


            print(f"{' '*7}With Error:\t\t{redText}{botchedSimulation['errorType']}{resetTextColor}")
            print(f"{' '*7}With Message:\t\t{redText}{botchedSimulation['errorMessage']}{resetTextColor}")
            print(f"{' '*7}At Line Number:\t\t{redText}{botchedSimulation['lineNumber']}{resetTextColor}")

            print(f"{tealColor}{' '*7}{'#'*4}{' '*7}Full Debug Traceback{' '*7}{'#'*4}{resetTextColor}")

            print(f"\t{orangeText}{'LINE NUMBER':<10}{yellowText}{'FUNCTION':>30}{resetTextColor}\t/path/to/crashed/{tealColor}script_name.py{resetTextColor}")
            print(f"\t{'---':<10}{'---':>30}\t{'---'}")
            for tracebackLine in botchedSimulation["fullTraceBack"]:
                scriptPath = tracebackLine.split(":")[0]
                scriptDir = p.dirname(scriptPath)
                scriptName = p.basename(scriptPath)
                lineNumber = tracebackLine.split(":")[1].split("in")[0].strip()
                functionName = tracebackLine.split(":")[1].split("in")[1].strip()
                print(f"\t{orangeText}{lineNumber:<10}{yellowText}{functionName:>30}{resetTextColor}\t{scriptDir}/{tealColor}{scriptName}{resetTextColor}")
    print(f"{redText}⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕")
    print(resetTextColor)
###########################################################################################

def print_prep_failed(errorMessage: str, stepName) -> None:
    redText = "\033[31m"
    yellowText = "\033[33m"
    resetTextColor = "\033[0m"


    wrappedErrorMessage = textwrap.fill(str(errorMessage), 80).replace("\n", "\n\t")

    print(redText+
          f"""
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕

    ██████╗ ██████╗ ███████╗██████╗  █████╗ ██████╗  █████╗ ████████╗██╗ ██████╗ ███╗   ██╗    
    ██╔══██╗██╔══██╗██╔════╝██╔══██╗██╔══██╗██╔══██╗██╔══██╗╚══██╔══╝██║██╔═══██╗████╗  ██║    
    ██████╔╝██████╔╝█████╗  ██████╔╝███████║██████╔╝███████║   ██║   ██║██║   ██║██╔██╗ ██║    
    ██╔═══╝ ██╔══██╗██╔══╝  ██╔═══╝ ██╔══██║██╔══██╗██╔══██║   ██║   ██║██║   ██║██║╚██╗██║    
    ██║     ██║  ██║███████╗██║     ██║  ██║██║  ██║██║  ██║   ██║   ██║╚██████╔╝██║ ╚████║    
    ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝   ╚═╝   ╚═╝ ╚═════╝ ╚═╝  ╚═══╝    

        ███████╗████████╗███████╗██████╗     ███████╗ █████╗ ██╗██╗     ███████╗██████╗        
        ██╔════╝╚══██╔══╝██╔════╝██╔══██╗    ██╔════╝██╔══██╗██║██║     ██╔════╝██╔══██╗       
        ███████╗   ██║   █████╗  ██████╔╝    █████╗  ███████║██║██║     █████╗  ██║  ██║       
        ╚════██║   ██║   ██╔══╝  ██╔═══╝     ██╔══╝  ██╔══██║██║██║     ██╔══╝  ██║  ██║       
        ███████║   ██║   ███████╗██║         ██║     ██║  ██║██║███████╗███████╗██████╔╝       
        ╚══════╝   ╚═╝   ╚══════╝╚═╝         ╚═╝     ╚═╝  ╚═╝╚═╝╚══════╝╚══════╝╚═════╝  
        At step {yellowText}{stepName}{redText}
        With Error:
        {yellowText}{wrappedErrorMessage}{redText}

        Please check the preparation log file for further details
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕{resetTextColor}
""")

    exit(1)

###########################################################################################


def print_performing_first_aid() -> None:
    
    yellowText = "\033[33m"
    resetTextColor = "\033[0m"

    print(yellowText+
          """
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕

       ██████  ███████ ██████  ███████  ██████  ██████  ███    ███ ██ ███    ██  ██████  
       ██   ██ ██      ██   ██ ██      ██    ██ ██   ██ ████  ████ ██ ████   ██ ██       
       ██████  █████   ██████  █████   ██    ██ ██████  ██ ████ ██ ██ ██ ██  ██ ██   ███ 
       ██      ██      ██   ██ ██      ██    ██ ██   ██ ██  ██  ██ ██ ██  ██ ██ ██    ██ 
       ██      ███████ ██   ██ ██       ██████  ██   ██ ██      ██ ██ ██   ████  ██████  
  
                 ███████ ██ ██████  ███████ ████████      █████  ██ ██████  
                 ██      ██ ██   ██ ██         ██        ██   ██ ██ ██   ██ 
                 █████   ██ ██████  ███████    ██        ███████ ██ ██   ██ 
                 ██      ██ ██   ██      ██    ██        ██   ██ ██ ██   ██ 
                 ██      ██ ██   ██ ███████    ██        ██   ██ ██ ██████  
                                                           
                             Your simulation has crashed...
                             drMD will attempt to rescue it...
         You may need to adjust your simulation parameters in your config file!
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕
          """
          +resetTextColor)
    

def print_first_aid_failed(errorOpenMM) -> None:
    """
    Prints the first aid failed message.

    Returns:    
        None
    """ 
    redText = "\033[31m"
    yellowText = "\033[33m"

    resetTextColor = "\033[0m"

    print(redText+
          f"""
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕

            ███████ ██ ██████  ███████ ████████      █████  ██ ██████  
            ██      ██ ██   ██ ██         ██        ██   ██ ██ ██   ██ 
            █████   ██ ██████  ███████    ██        ███████ ██ ██   ██ 
            ██      ██ ██   ██      ██    ██        ██   ██ ██ ██   ██ 
            ██      ██ ██   ██ ███████    ██        ██   ██ ██ ██████  
                                                                    
                                                                    
                ███████  █████  ██ ██      ███████ ██████              
                ██      ██   ██ ██ ██      ██      ██   ██             
                █████   ███████ ██ ██      █████   ██   ██             
                ██      ██   ██ ██ ██      ██      ██   ██             
                ██      ██   ██ ██ ███████ ███████ ██████

        drMD failed to rescue your simulation.
        The following error was returned by OpenMM:

        {yellowText}{errorOpenMM}{redText}

        This is likely due to unphysically high energies in your simulation.
    
        Try some of the following:
        -> reduce the timestep of your simulation
        -> reduce the temperature of your simulation
        -> if you are using any restraints, reduce their force constants

⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕
                                                          
"""
+resetTextColor)
    


###########################################################################################
def print_config_error(configDisorders) -> None:
    """
    Prints an error message indicating that the config file was not found.

    Returns:
        None
    """
    redText = "\033[31m"
    yellowText = "\033[33m"
    orangeText = "\033[38;5;208m"
    greenText = "\033[32m"
    resetTextColor = "\033[0m"


    print(redText+
          """
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕
          
 ▄████▄  ▒█████   ███▄    █   █████▒██▓ ▄████    ▓█████  ██▀███   ██▀███   ▒█████   ██▀███  
▒██▀ ▀█ ▒██▒  ██▒ ██ ▀█   █ ▓██   ▒▓██▒██▒ ▀█▒   ▓█   ▀ ▓██ ▒ ██▒▓██ ▒ ██▒▒██▒  ██▒▓██ ▒ ██▒
▒▓█    ▄▒██░  ██▒▓██  ▀█ ██▒▒████ ░▒██▒██░▄▄▄░   ▒███   ▓██ ░▄█ ▒▓██ ░▄█ ▒▒██░  ██▒▓██ ░▄█ ▒
▒▓▓▄ ▄██▒██   ██░▓██▒  ▐▌██▒░▓█▒  ░░██░▓█  ██▓   ▒▓█  ▄ ▒██▀▀█▄  ▒██▀▀█▄  ▒██   ██░▒██▀▀█▄  
▒ ▓███▀ ░ ████▓▒░▒██░   ▓██░░▒█░   ░██░▒▓███▀▒   ░▒████▒░██▓ ▒██▒░██▓ ▒██▒░ ████▓▒░░██▓ ▒██▒
░ ░▒ ▒  ░ ▒░▒░▒░ ░ ▒░   ▒ ▒  ▒ ░   ░▓  ░▒   ▒    ░░ ▒░ ░░ ▒▓ ░▒▓░░ ▒▓ ░▒▓░░ ▒░▒░▒░ ░ ▒▓ ░▒▓░
  ░  ▒    ░ ▒ ▒░ ░ ░░   ░ ▒░ ░      ▒ ░ ░   ░     ░ ░  ░  ░▒ ░ ▒░  ░▒ ░ ▒░  ░ ▒ ▒░   ░▒ ░ ▒░
░       ░ ░ ░ ▒     ░   ░ ░  ░ ░    ▒ ░ ░   ░       ░     ░░   ░   ░░   ░ ░ ░ ░ ▒    ░░   ░ 
░ ░         ░ ░           ░         ░       ░       ░  ░   ░        ░         ░ ░     ░     
░          
          """)
    print(f"{resetTextColor}The following disorders have been found in your config file:")    
    print(f"{resetTextColor}Colour Key: | {greenText}Input Correct{resetTextColor} | {orangeText}Non Fatal, Default Used{resetTextColor} | {redText}Fatal Issue{resetTextColor} |")    


    def print_config_text(argName, argDisorder, textColor, indentationLevel=0) -> None:
        print(f"{' '*(indentationLevel*3+2)}{yellowText}{argName}: {textColor}{argDisorder}{resetTextColor}")
    
    def loop_disorder_dict(argName, disorderDict, indentationLevel=0) -> None:
        print(f"{'--'*indentationLevel}--> In sub-entry {yellowText}{argName}{resetTextColor}:")
        for argName, argDisorder in disorderDict.items():
            if argDisorder is None:
                print_config_text(argName, argDisorder, greenText, indentationLevel)
            elif isinstance(argDisorder, str):
                if "default" in argDisorder.lower():
                    print_config_text(argName, argDisorder, orangeText, indentationLevel)
                else:
                    print_config_text(argName, argDisorder, redText, indentationLevel)
            elif isinstance(argDisorder, list):
                print_config_text(argName, argDisorder, redText, indentationLevel)
            elif isinstance(argDisorder, dict):
                loop_disorder_dict(argName, argDisorder, indentationLevel + 1)


    for infoName, infoDisorders in configDisorders.items():
        if infoDisorders is None:
            continue
        print(f"> For the config entry {yellowText}{infoName}{resetTextColor}, the following problems were found:")
        for argName, argDisorder in infoDisorders.items():
            if argDisorder is None:
                print_config_text(argName, argDisorder, greenText, 0 )
            elif isinstance(argDisorder, str):
                if "default" in argDisorder.lower():
                    print_config_text(argName, argDisorder, orangeText, 0 )
                else:
                    print_config_text(argName, argDisorder, redText, 0)
            elif isinstance(argDisorder, list):
                    print_config_text(argName, argDisorder, redText, 0)
            elif isinstance(argDisorder, dict):
                loop_disorder_dict(argName, argDisorder)

    print(f"{redText}⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕")

    print(resetTextColor)
    exit(1)




###########################################################################################
def print_pdb_error(pdbDisorders) -> None:
    """
    Prints an error message indicating that problems have been found in the pdb files

    Returns:
        None
    """
    pdbFixerLink = "https://github.com/openmm/pdbfixer"
    alphafold3Link = "https://alphafoldserver.com/"

    ionResNames = list(drListInitiator.get_ion_residue_names())

    ionResNameText = ""
    for i in range(0, len(ionResNames), 10):
        ionResNameText += '\t\t' + ', '.join(ionResNames[i:i+10]) + '\n'

    pdbRules = {
        "01_broken_protein_chains" : 
    {
    "number": "1",
    "text": "Protein chains must be contiguous (i.e., must not have missing residues)",
    "help": f"To fix this, drMD recommends you rebuild your PDB using OpenMM's pdbFixer: {pdbFixerLink}"
    },
    "02_residues_missing_atoms" :
    {
    "number": "2",
    "text": "All amino acid residues must have the expected number of sidechain atoms",
    "help": f"""If you have a crystal structure, this may be the result of low electron density
\t\tfor some residue sidechains. To fix this, try getting a structural prediction of your protein from 
\t\tAlphafold's webserver {alphafold3Link}"""
    },
    "03_residues_with_duplicate_atoms" :
    {
    "number": "3",
    "text": "Residues must not have duplicate atoms",
    "help": """If you have a crystal structure, you may have multiple conformers of some sidechains,
\t\tTo fix this, use the following pymol commands:\n\t\tremove not alt \'\'+A\n\t\talter all, alt = \'\' """
    },
    "04_atoms_with_no_chain_id" :
    {
    "number": "4",
    "text": "All atoms must have a chain identifier",
    "help": "This can be done with ease in your text editor, or with pdbUtils"
    },
    "05_protein_chains_unique_chain_ids" :
    {
    "number": "5",
    "text": "Covalently bound protein chains must each have a unique chain identifier",
    "help": "Change the chain identifier for the protein chains"
    },
    "06_ligands_and_protein_sharing_chain_ids" :
    {
    "number": "6",
    "text": "Protein chains cannot contain ligand residues",
    "help": "Change the chain identifier for the ligand residues"
    },
    "07_organometallic_ligands" :
    {
    "number": "7",
    "text": "No organometallic residues",
    "help": """Currently, drMD does not support organometallic ligands,
\t\tyou may be able to separate your ligand into the organic moiety and ion components"""
    },
    "08_non-canonical_amino_acids" :
    {
    "number": "8",
    "text": "No non-canonical amino acid residues",
    "help": """Currently, drMD does not support non-canonical amino acid residues. 
\t\tWe are currently working on a solution to this"""
    },
    "09_ions_with_incorrect_names" :
    {
    "number": "9",
    "text": "Ions must have the residue and atom names compatible with AMBER",
    "help": f"""Set your ion residue and atom name to be one of these:\n {ionResNameText}"""
    }
    }


    tealText = "\033[38;5;37m" 
    redText = "\033[31m"
    yellowText = "\033[33m"
    orangeText = "\033[38;5;208m"
    greenText = "\033[32m"
    resetTextColor = "\033[0m"


    print(redText+
          """
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕

            ██▓███  ▓█████▄  ▄▄▄▄      ▓█████  ██▀███   ██▀███   ▒█████   ██▀███  
            ▓██░  ██▒▒██▀ ██▌▓█████▄    ▓█   ▀ ▓██ ▒ ██▒▓██ ▒ ██▒▒██▒  ██▒▓██ ▒ ██▒
            ▓██░ ██▓▒░██   █▌▒██▒ ▄██   ▒███   ▓██ ░▄█ ▒▓██ ░▄█ ▒▒██░  ██▒▓██ ░▄█ ▒
            ▒██▄█▓▒ ▒░▓█▄   ▌▒██░█▀     ▒▓█  ▄ ▒██▀▀█▄  ▒██▀▀█▄  ▒██   ██░▒██▀▀█▄  
            ▒██▒ ░  ░░▒████▓ ░▓█  ▀█▓   ░▒████▒░██▓ ▒██▒░██▓ ▒██▒░ ████▓▒░░██▓ ▒██▒
            ▒▓▒░ ░  ░ ▒▒▓  ▒ ░▒▓███▀▒   ░░ ▒░ ░░ ▒▓ ░▒▓░░ ▒▓ ░▒▓░░ ▒░▒░▒░ ░ ▒▓ ░▒▓░
            ░▒ ░      ░ ▒  ▒ ▒░▒   ░     ░ ░  ░  ░▒ ░ ▒░  ░▒ ░ ▒░  ░ ▒ ▒░   ░▒ ░ ▒░
            ░░        ░ ░  ░  ░    ░       ░     ░░   ░   ░░   ░ ░ ░ ░ ▒    ░░   ░ 
                        ░     ░            ░  ░   ░        ░         ░ ░     ░  
        drMD has found problems with your PDB files that will cause prep steps to fail
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕
          """
          +resetTextColor)


    print(f"{tealText}\tdrMD requires PDB files to follow these rules:\n{resetTextColor}")

    for disorderTag, problemPdbs in pdbDisorders.items():
        disorder = pdbRules[disorderTag]
        print(f"\t{tealText}{disorder['number']}{resetTextColor}.{' '*2}{disorder['text']}")
        if len(problemPdbs) == 0:
            print(f"\t\t{greenText}Check passed!{resetTextColor}")
            continue
        if len(problemPdbs) < 6:
            print(f"{yellowText}\t\tThe following PDB files failed this check:{resetTextColor}")
            print(f"{redText}\t\t{' '.join(problemPdbs)}")
        else:
            print(f"{redText}\t\t{len(problemPdbs)}{yellowText} failed this check")
        print(f"\t\t{greenText}{disorder['help']}")

    print(f"""{redText}
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕
    {resetTextColor}""")
    exit(1)

###########################################################################################



if __name__ == "__main__":
    print_drMD_logo()
    print_pdb_error()
    print_performing_first_aid()

    print_config_error({
        "pathInfo": {
            "inputDir": "inputDir must be a path to a directory containing .pdb files",
            "outputDir": None
            },
        "hardwareInfo": {
            "parallelCPU": "parallelCPU must be an int greater than 1",
            "platform": "platform must be CUDA, OpenCL, or CPU",
            "subprocessCpus": "no entry found, using a default of 1"
        }
    })
  