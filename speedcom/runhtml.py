from http.server import BaseHTTPRequestHandler, HTTPServer
import urllib
import json
from os import path
from urllib.parse import urlparse
#import function_name as pyf


curdir = './' # current directory
sep = '/'
# MIME-TYPE
mimedic = [
    ('.html', 'text/html'),
    ('.htm', 'text/html'),
    ('.js', 'application/javascript'),
    ('.css', 'text/css'),
    ('.json', 'application/json'),
    ('.png', 'image/png'),
    ('.jpg', 'image/jpeg'),
    ('.gif', 'image/gif'),
    ('.txt', 'text/plain'),
    ('.avi', 'video/x-msvideo'),
]


class ServerHTTP(BaseHTTPRequestHandler):
    # http get request
    def do_GET(self):
        sendReply = False
        querypath = urlparse(self.path)
        filepath, query = querypath.path, querypath.query

        if filepath.endswith('/'):
            filepath += 'speedcom.html'
        filename, fileext = path.splitext(filepath)
        for e in mimedic:
            if e[0] == fileext:
                mimetype = e[1]
                sendReply = True

        if sendReply == True:
            try:
                with open(path.realpath(curdir + sep + filepath),'rb') as f:
                    content = f.read()
                    self.send_response(200)
                    self.send_header('Content-type',mimetype)
                    self.end_headers()
                    self.wfile.write(content)
            except IOError:
                self.send_error('File Not Found: ')


    # http post request
    def do_POST(self):
        path = self.path
        # get request path. -> http://domain:port{path} <-this path
        print(path)

        # fetch post body data
        # we use application/json protocol
        mpath,margs=urllib.parse.splitquery(self.path)
        source = self.rfile.read(int(self.headers['content-length']))
        print(source)
        data = json.loads(source)


        ######
        #
        # Run the python function
        #
        #
        ######
        if path == '/input':
            print(data["input"])
            # the name of the python function
            
            # read the value of each characteristics
            fvalue = open("output/characteristics.txt", 'r')
            lvalue = fvalue.readlines()
            fvalue.close()
            lvalue = lvalue[-1].split('\t')
            print(lvalue)
            value = []
            for i in range(len(lvalue)):
                if lvalue[i] == '\n':
                    pass
                else:
                    value.append(eval(lvalue[i]))
            print(value)


        # response data to html client.
        self.send_response(200)
        self.send_header("Content-type","Application/json")
        self.send_header("Access-Control-Allow-Origin", "*");
        self.end_headers()


        # use json.dumps to format
        response = json.dumps({
            "input" : data["input"],
            "value" : value
        }).encode()
        # send
        self.wfile.write(response)


def start_server(port):
    print("listened at ", port)
    http_server = HTTPServer(('', int(port)), ServerHTTP)
    http_server.serve_forever()


if __name__ == "__main__":
    start_server(8000)
