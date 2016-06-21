import pickle
from os.path import join, exists

CACHE_DIR = 'cache'

SERIALIZER = pickle.dumps
DESERIALIZER = pickle.loads

def cache_file_for_(fct, args, kwargs):
    return join(
        CACHE_DIR,
        '{fct}_{args}_{kwargs}'.format(
            fct=fct.__name__,
            args=','.join(map(str, args)),
            kwargs=','.join(['{0}={1}'.format(key, value) for (key, value) in sorted(kwargs.items())]),
        ),
    )

def cached(fct, args, kwargs):
    cache_file_path = cache_file_for_(fct, args, kwargs)
    if exists(cache_file_path):
        return DESERIALIZER(open(cache_file_path).read())
    else:
        fct_return = fct(*args, **kwargs)
        with open(cache_file_path, 'w') as fh:
            fh.write(SERIALIZER(fct_return))
        return fct_return

if __name__ == '__main__':
    def add(a, b, debug=False):
        return a + b

    print cached(add, (1, 2), dict(debug=True))
    print cached(add, (1, 2), dict(debug=True))
