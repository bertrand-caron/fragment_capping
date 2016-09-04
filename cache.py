import pickle
from os.path import join, exists
from hashlib import md5

CACHE_DIR = 'cache'

SERIALIZER = pickle.dumps
DESERIALIZER = pickle.loads

def signature(fct_args, fct_kwargs, hashed):
    if hashed:
        signature_fct = lambda x: md5(x.encode()).hexdigest()
    else:
        signature_fct = lambda x: x

    return signature_fct(
        ','.join(map(str, fct_args))
        +
        '_'
        +
        ','.join(['{0}={1}'.format(key, value) for (key, value) in sorted(fct_kwargs.items())])
        +
        '.pkl'
    )

def cache_file_for_(fct, fct_args, fct_kwargs, hashed=False):
    return join(
        CACHE_DIR,
        '{fct}_{signature}'.format(
            fct=fct.__name__,
            signature=signature(fct_args, fct_kwargs, hashed),
        ),


    )

def cached(fct, fct_args, fct_kwargs, hashed=False):
    cache_file_path = cache_file_for_(fct, fct_args, fct_kwargs, hashed=hashed)
    if exists(cache_file_path):
        return DESERIALIZER(open(cache_file_path, 'rb').read())
    else:
        fct_return = fct(*fct_args, **fct_kwargs)
        with open(cache_file_path, 'wb') as fh:
            fh.write(SERIALIZER(fct_return))
        return fct_return

if __name__ == '__main__':
    def add(a, b, debug=False):
        return a + b

    print(cached(add, (1, 2), dict(debug=True)))
    print(cached(add, (1, 2), dict(debug=True)))
    print(cached(add, (1, 2), dict(debug=True), hashed=True))
    print(cached(add, (1, 2), dict(debug=True), hashed=True))
