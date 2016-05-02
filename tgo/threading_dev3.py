from multiprocessing import Pool

def f(x):
    return x*x

if __name__ == '__main__':
    #p = Pool(5)
    p = Pool()
    #Pool(
    print(p.map(f, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                    11, 12, 13, 14, 15, 16, 17, 980e3, 1, 2]))

# from multiprocessing import Process, Queue
#
#
# def f(q):
#     q.put([42, None, 'hello'])
#
#
# if __name__ == '__main__':
#     q = Queue()
#     p = Process(target=f, args=(q,))
#     p.start()
#     print(q.get() ) # prints "[42, None, 'hello']"
#     p.join()
#
#
# from multiprocessing import Process, Queue
#
# def f(q):
#     q.put([42, None, 'hello'])
#
# if __name__ == '__main__':
#     q = Queue()
#     p = Process(target=f, args=(q,))
#     p.start()
#     print(q.get())    # prints "[42, None, 'hello']"
#     p.join()