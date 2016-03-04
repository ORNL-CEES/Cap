from .PyCap import PropertyTree
from weakref import ref

__all__ = ['Observer', 'Observable', 'Experiment']

class Observer(object):
    '''Observer.
    
    An observer provides a representation of an observable, called the subject.
    It may be attached to one or several subjects.
    Subjects are able to notify the observers that have been attached to them
    when some change has occured and an update of their representation is
    necessary.
    
    Examples
    --------
    >>> ptree = PropertyTree()
    >>> ptree.put_string('type', 'ConcreteObserver')
    >>> observer = Observer(ptree)

    >>> observer.update(subject)

    >>> subject.attach(observer)
    >>> subject.notify()

    See Also
    --------
    Observable
    '''

    _builders = {}

    def __new__(cls, ptree, *args, **kwargs):
        '''Create a new instance of class cls.
        
        Parameters
        ----------
        cls : class
        ptree : PropertyTree
        
        Raises
        ------
        KeyError
            The observer type is not valid.
        RuntimeError
            The property tree is missing the key 'type' or is not able to
            convert its value to string.
        '''
        t = ptree.get_string('type')
        return Observer._builders[t].__new__(cls, ptree, *args,  **kwargs)

    def update(self, subject, *args, **kwargs):
        '''Update the representation of the Observable.
        
        Must be overriden in derived classes.
        
        Parameters
        ----------
        subject : Observable
        args : list
            Arbitrary argument list.
        kwargs : dict
            Keyword arguments.
        '''
        raise RuntimeError('Method Observer.update(...) must be overloaded')



class Observable(object):
    '''Observable.
    
    Maintains a list of observers and notifies them when something has changed.
    
    Attributes
    ----------
    _observers : list
        Stores weak references to the observers.

    Examples
    --------
    >>> ptree = PropertyTree()
    >>> ptree.put_string('type', 'ConcreteObservable')
    >>> subject = Observable(ptree)

    >>> observer.update(subject)

    >>> subject.attach(observer)
    >>> subject.notify()

    See Also
    --------
    Observer
    '''

    _builders = {}

    def __new__(cls, ptree, *args, **kwargs):
        '''Create a new instance of class cls.
        
        Parameters
        ----------
        cls : class
        ptree : PropertyTree
        
        Raises
        ------
        KeyError
            The observer type is not valid.
        RuntimeError
            The property tree is missing the key 'type' or is not able to
            convert its value to string.
        '''
        t = ptree.get_string('type')
        return Observable._builders[t].__new__(cls, ptree, *args, **kwargs)
    def __init__(self):
        self._observers = []

    def attach(self, observer):
        '''Attach an observer.
        
        Adds a weak reference to the observer to the list of observers.
        
        Parameters
        ----------
        observer : Observer
        
        Raises
        ------
        RuntimeError
            You are trying to attach the same observer twice.
        '''
        for weak_reference in self._observers:
            if weak_reference() is observer:
                raise RuntimeError('Observer already registered')
        self._observers.append(ref(observer))

    def detach(self, observer):
        '''Detach an observer.
        
        Removes the observer from the list.
        
        Parameters
        ----------
        observer : Observer
        
        Raises
        ------
        RuntimeError
            You are trying to detach an observer that is not attached.
        '''
        for weak_reference in self._observers:
            if weak_reference() is observer:
                self._observers.remove(weak_reference)
                return
        raise RuntimeError('Observer not attached')
        
    def notify(self, *args, **kwargs):
        '''Notify all registered observers that some change has occured.
        '''
        for weak_reference in self._observers:
            observer = weak_reference()
            if observer is not None:
                observer.update(self, *args, **kwargs)
            else:
                raise RuntimeError("Observer is no longer alive")


class Experiment(Observable):
    '''Experiment

    Base class for the electrochemical measurement techniques in PyCap.

    Attributes
    ----------
    _data : dict
        Holds the measurement information.
    '''
    def __init__(self):
        Observable.__init__(self)
        self._data = {}
    def run(self, device):
        '''Run the experiment on an energy storage device.

        Parameters
        ----------
        device : EnergyStorageDevice
        '''
        raise RuntimeError('Method Experiment.run(...) must be overloaded')
    
