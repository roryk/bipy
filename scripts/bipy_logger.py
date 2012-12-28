import zmq
import sys

if __name__ == "__main__":
    connect_to = sys.argv[1]
    log_file = sys.argv[2]

    out_handle = open(log_file, "w")
    out_handle.write("Starting logging.")
    out_handle.flush()

    context = zmq.Context()
    b = context.socket(zmq.SUB)
    b.setsockopt(zmq.SUBSCRIBE, '')
    b.connect(connect_to)

    while True:
        message = b.recv()
        out_handle.write(message + "\n")
        out_handle.flush()
