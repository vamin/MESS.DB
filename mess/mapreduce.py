# -*- coding: utf-8 -*-
# Copyright 2013-2014 Victor Amin, http://vamin.net/

"""MESS.DB mapreduce module

This module implements the mapreduce paradigm by extending classes from the
mincemeatpy module.
"""

import os
import socket
import sys
import time
from collections import OrderedDict

sys.path.insert(1, os.path.join(os.path.dirname(__file__), 'mincemeatpy'))
import mincemeat

from mess.log import Log

# override defualt logging use with custom Log class
logging = Log('console')
mincemeat.logging = Log('console')

DEFAULT_PORT = mincemeat.DEFAULT_PORT


class Client(mincemeat.Client, object):
    """Extends mincemeat.Client to support 'declined' command"""
    
    def __init__(self):
        """Call parent init."""
        super(Client, self).__init__()
    
    def handle_connect(self):
        """Keep handle_connect() pass."""
        pass
    
    def call_mapfn(self, command, data):
        """Call mapfn, or decline if mapfn not set."""
        if not self.mapfn:
            time.sleep(1)
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
        """Call reducefn, or decline and close if reducefn not set."""
        if not self.reducefn:
            time.sleep(1)
            self.send_command('declined')
            self.handle_close()
        else:
            logging.info("Reducing %s" % str(data[0]))
            results = self.reducefn(data[0], data[1])
            self.send_command('reducedone', (data[0], results))
    
    def run(self, hostname, port, sleep=2):
        """Run client."""
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
    """Extends mincemeat.Server to add 'run' method and to use extended
    ServerChannel class"""
    
    def __init__(self):
        """Call parent init."""
        super(Server, self).__init__()
    
    def run(self, password=None, port=DEFAULT_PORT, debug=0):
        """Run server."""
        if password is None:
            password = self.password
        return self.run_server(password, port)
    
    def handle_accept(self):
        """Accept a connection, uses extended ServerChannel."""
        conn, addr = self.accept()
        sc = ServerChannel(conn, self.socket_map, self)
        sc.password = self.password


class ServerChannel(mincemeat.ServerChannel):
    """Extends mincemeat.ServerChannel to to support 'declined' command"""
    
    def __init__(self, conn, map, server):
        """Reimplement parent init."""
        mincemeat.Protocol.__init__(self, conn, map=map)
        self.server = server
        self.start_auth()
    
    def declined(self, command, data):
        """Go to next task if a task is declined."""
        self.start_new_task()
    
    def process_command(self, command, data=None):
        """Process 'mapdone', 'reducedone', and 'declined' commands."""
        commands = {
            'mapdone': self.map_done,
            'reducedone': self.reduce_done,
            'declined': self.declined
            }

        if command in commands:
            commands[command](command, data)
        else:
            mincemeat.Protocol.process_command(self, command, data)
