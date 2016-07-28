from .genapi import GenAPI


def connect(username, password, url, server_type):
    """ Connect to Resolwe server
    Args:
        username (:obj:`str`):
        password (:obj:`str`):
        url (:obj:`str`): url of the server you are connecting
        server_type (:obj:`str`): genesis or resolwe

    Returns:
        Instance of GenAPI or ResolweAPI

    """

    if server_type == 'genesis':
        try:
            return GenAPI(username, password, url)
        except Exception as e:
            raise ResolweAuthException(e.args[0]) from e
    else:
        """ Not yet supported """
        pass


class ResolweAuthException(Exception):
    """A login error occurred."""



