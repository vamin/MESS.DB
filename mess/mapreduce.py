import logging
import socket
import sys
import time
from collections import OrderedDict

from mess.mincemeatpy import mincemeat

DEFAULT_PORT = mincemeat.DEFAULT_PORT

class Client(mincemeat.Client):
    def __init__(self):
        super(Client, self).__init__()
    
    def handle_error(self):
        raise
    
    def handle_connect(self):
        pass
    
    def call_mapfn(self, command, data):
        if not self.mapfn:
            self.send_command('declined')
        else:
            logging.info("Mapping %s" % str(data[0]))
            results = OrderedDict()
            for k, v in self.mapfn(data[0], data[1]):
                if k not in results:
                    results[k] = []
                results[k].append(v)
            if self.collectfn:
                for k in results:
                    results[k] = [self.collectfn(k, results[k])]
            self.send_command('mapdone', (data[0], results))
    
    def call_reducefn(self, command, data):
        if not self.reducefn:
            self.send_command('declined')
        else:
            logging.info("Reducing %s" % str(data[0]))
            results = self.reducefn(data[0], data[1])
            self.send_command('reducedone', (data[0], results))
    
    def run(self, hostname, port, sleep=2, debug=0):  
        if debug:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)
        while True:
            try:
                self.conn(hostname, port)
                break
            except socket.error:
                exc_info = sys.exc_info()
                logging.debug('%s.{hostname=%s, port=%s}:%s',
                              exc_info[0],
                              hostname,
                              port,
                              exc_info[1])
                time.sleep(sleep)
                self.__init__()
            except KeyboardInterrupt:
                return
            except:
                exc_info = sys.exc_info()
                logging.exception('%s:%s', exc_info[0], exc_info[1])
                break


class Server(mincemeat.Server):
    def __init__(self):
        super(Server, self).__init__()
    
    def run(self, password=None, port=DEFAULT_PORT, debug=0):
        if password is None:
            password = self.password
        if debug:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)
        self.run_server(password, port)


class ServerChannel(mincemeat.ServerChannel):
    def __init__(self):
        super(ServerChannel, self).__init__()
    
    def declined(self, command, data):
        self.start_new_task()
    
    def process_command(self, command, data=None):
        commands = {
            'mapdone': self.map_done,
            'reducedone': self.reduce_done,
            'declined': self.declined
            }

        if command in commands:
            commands[command](command, data)
        else:
            Protocol.process_command(self, command, data)
