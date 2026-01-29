from .rips_persistence import RipsPersistence

__all__ = [
    'RipsPersistence',
]

try:
    # If no array_api_compat
    from .cubical_persistence import CubicalPersistence
    __all__ += [
        'CubicalPersistence',
    ]
except ModuleNotFoundError:
    pass

try:
    # if no CGAL
    from .cech_persistence import CechPersistence, WeightedCechPersistence
    
    __all__ += [
        'CechPersistence',
        'WeightedCechPersistence',
    ]
except ImportError:
    pass
