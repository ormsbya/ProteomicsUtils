import logging
import sys, os

def logger_config(logger_name, logPath=False, file_log=True, print_log=True):
    """
    Takes logging.info and above messages, printing them to console and saving to a text file

    Parameters:
    logger_name: string
        name of the logger to be instantiated e.g. analysis
    logPath: string or False
        file path to save the Log file to, defaults to false in which case Log file is saved to current working directory
    file_log: bool
        default True, instantiates the creation of the logger file. If false, does not save log to file.
    print_log: bool
        default True, instantiates the messages to be printed to console. If false, does not print.

    """
    #to set output path as current working directory if LogPath not given
    if not logPath:
        logPath = os.getcwd()
    else:
        pass
    #set file name, create formatter for log output
    fileName = 'Log file'
    logFormatter = logging.Formatter("%(asctime)s %(name)s: [%(levelname)-5.5s]  %(message)s")

    currentLogger = logging.getLogger(logger_name)
    currentLogger.setLevel(logging.DEBUG)
    #check if current logger already has handlers attached
    if len(currentLogger.handlers):
        return currentLogger
    else:
        if file_log:
            fileHandler = logging.FileHandler("{0}/{1}.log".format(logPath, fileName))
            fileHandler.setFormatter(logFormatter)
            currentLogger.addHandler(fileHandler)
        if print_log:
            consoleHandler = logging.StreamHandler()
            consoleHandler.setFormatter(logFormatter)
            currentLogger.addHandler(consoleHandler)

        return currentLogger
