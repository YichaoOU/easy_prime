"""
Using __init__.py to organize the structure.
"""
try:
	from ._target_mutation import target_mutation
except ImportError as e:
	print (e)

__all__ = ["target_mutation"]
__version__ = "1.1"