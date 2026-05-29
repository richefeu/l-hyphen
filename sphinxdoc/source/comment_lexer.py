from pygments.lexer import RegexLexer
from pygments.token import Comment

class CommentLexer(RegexLexer):
    """
    Lexer for highlighting comments in a text file.
    """
    name = 'CommentLexer'
    aliases = ['commentlexer']

    tokens = {
        'root': [
            (r'#.*$', Comment.Single),  # Matches comments starting with #
            (r'!.*$', Comment.Single),  # Matches comments starting with !
        ]
    }
