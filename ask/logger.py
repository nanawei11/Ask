import logging

class Logger:
    def __init__(self):
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)
        self.logger.handlers.clear()
        fmt = logging.Formatter('[%(asctime)s] [%(levelname)s] %(message)s', '%Y-%m-%d %H:%M:%S')
        logg = logging.StreamHandler()
        logg.setFormatter(fmt)
        self.logger.addHandler(logg)

    def info(self, message):
        self.logger.info(message)

    def war(self, message):
        self.logger.warn(message)

    def error(self, message):
        self.logger.error(message)