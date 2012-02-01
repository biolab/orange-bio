"""
Caching framework for cached kegg api calls.
 
"""
import os
import UserDict
import sqlite3
import cPickle as pickle

from datetime import datetime


class Store(object):
    def __init__(self):
        self.timestamp = 0
        
    def open(self):
        raise NotImplementedError
    
    def __enter__(self):
        return self
    
    def __exit__(self, *args):
        pass
    
    
class ShelveStore(Store):
    def __init__(self):
        import shelve
        self.fname = shelve.open


class Sqlite3Store(Store, UserDict.DictMixin):
    def __init__(self, filename):
        self.filename = filename
        self.con = sqlite3.connect(filename)
        self.con.execute("""
        CREATE TABLE IF NOT EXISTS cache 
            (key TEXT UNIQUE,
             value TEXT
            )
        """)
        self.con.commit()
        
    def __getitem__(self, key):
        cur = self.con.execute("""
            SELECT value
            FROM cache
            WHERE key=?
        """, (key,))
        r = cur.fetchall()
        
        if not r:
            raise KeyError(key)
        else:
            return pickle.loads(str(r[0][0]))
    
    def __setitem__(self, key, value):
        value = pickle.dumps(value)
        self.con.execute("""
            INSERT INTO cache
            VALUES (?, ?)
        """, (key, value))
        self.con.commit()
        
    def __delitem__(self, key):
        self.con.execute("""
            DELETE FROM cache
            WHERE key=?
        """, (key,))
        self.con.commit()
        
    def keys(self):
        cur = self.con.execute("""
            SELECT key
            FROM cache
        """)
        return [str(r[0]) for r in cur.fetchall()]
        
    def close(self):
        pass
    
    
class DictStore(Store, UserDict.DictMixin):
    def __init__(self):
        Store.__init__(self)
        
    def close(self):
        pass
    
    
from functools import wraps
from contextlib import closing


class cache_entry(object):
    def __init__(self, value, mtime=None, expires=None):
        self.value = value
        self.mtime = mtime
        self.expires = expires
        

class cached_wrapper(object):
    """ TODO: needs documentation
    """
    def __init__(self, function, instance, class_, cache_store, last_modified=None):
        self.function = function
        self.instance = instance
        self.class_ = class_
        self.cache_store = cache_store
        self.last_modified = last_modified
        
    def has_key(self, key):
        with closing(self.cache_store()) as store:
            return key in store
    
    def key_from_args(self, args, kwargs=None):
        key = self.function.__name__ + repr(args)
        return key
    
    def invalidate_key(self, key):
        with closing(self.cache_store()) as store:
            del store[key]
            
    def last_modified_from_args(self, args, kwargs=None):
        key = self.key_from_args(args, kwargs)
        if self.instance is not None:
            self.instance.last_modified(args)
        
    def invalidate_args(self, args):
        return self.invalidate_key(self.key_from_args(args))
        
    def invalidate_all(self):
        prefix = self.key_from_args(()).rstrip(",)")
        with self.cache_store() as store:
            for key in store.keys():
                if key.startswith(prefix):
                    del store[key]
    
    def memoize(self, args, kwargs, value, timestamp=None):
        key = self.key_from_args(args, kwargs)
        if timestamp is None:
            timestamp = datetime.now()
            
        with closing(self.cache_store()) as store:
            store[key] = cache_entry(value, mtime=timestamp)
        
    def __call__(self, *args):
        key = self.key_from_args(args)
        with closing(self.cache_store()) as store:
            valid = True
            if key not in store:
                valid = False
            else:
                entry = store[key]
                rval = entry.value
                
                if not self.is_entry_valid(entry, args):
                    valid = False
                    
            if not valid:
                rval = self.function(self.instance, *args)
                store[key] = cache_entry(rval, datetime.now(), None)
        
        return rval
        
    def min_timestamp(self, args):
        key = self.key_from_args(args)
        return datetime.fromtimestamp(0)
    
    def is_entry_valid(self, entry, args):
        return self.min_timestamp(args) < entry.mtime 
    
        
class cached_method(object):
    def __init__(self, function):
        self.function = function
        
    def __get__(self, instance, owner):
        if instance is not None:
            return cached_wrapper(self.function, instance, owner,
                                  self.get_cache_store(instance, owner))
        return self
    
    def get_cache_store(self, instance, owner):
        if hasattr(instance, "cache_store"):
            return instance.cache_store
        elif not hasattr("_cached_method_cache"):
            instance._cached_method_cache = DictStore()
        return instance._cached_method_cache


class bget_cached_method(cached_method):
    def __get__(self, instance, owner):
        if instance is not None:
            return cached_wrapper(self.function, instance, owner,
                                  self.get_cache_store(instance, owner),
                                  self.get_last_modified(instance, owner))
        return self
    
    def get_last_modified(self, instance, owner):
        if hasattr(instance, "last_modified"):
            return instance.last_modified
    
def touch_dir(path):
    path = os.path.expanduser(path)
    if not os.path.exists(path):
        os.makedirs(path)
    