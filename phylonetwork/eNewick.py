from pyparsing import Combine, Optional, Literal, CaselessLiteral, \
    Word, alphanums, alphas, alphas8bit, \
    nums, oneOf, Group, Dict, Forward, \
    ParseResults, CharsNotIn, ZeroOrMore

cvtInt = lambda s,l,toks: int(toks[0])
cvtReal = lambda s,l,toks: float(toks[0])

def makeeNewickParser():
    # atoms    
    lparen    = Literal("(").suppress()
    rparen    = Literal(")").suppress()
    colon     = Literal(":").suppress()
    semicolon = Literal(";").suppress()
    comma     = Literal(",").suppress()
    point     = Literal(".")
    e         = CaselessLiteral("E")
    sharp     = Literal("#").suppress()

    # terminal
    name    = Word(alphanums + alphas8bit + "_" + "-" + "." + "+" + "&" + "/" + "~" + "{" + "}" + "*" + "'" + '"' + '\\' + '?')
    string  = Word(alphas)
    fnumber = Combine(
        Word("+-"+nums, nums) + 
        Optional(point + Optional(Word(nums))) +
        Optional(e + Word("+-"+nums, nums))
        ).setParseAction(cvtReal)
    number = Combine(
        Word(nums)).setParseAction(cvtInt)

    label = \
        Optional(name).setResultsName("label") + \
        Optional(
            sharp +
            Optional(string).setResultsName("type") +
            number.setResultsName("tag")
            ) + \
        Optional(colon + fnumber).setResultsName("length")

    subtree = Forward()
    subtreelist = Forward()
    subtree << \
    Group(((lparen + subtreelist + rparen).setResultsName("subtree") |
           label
           ) + Optional(label)
          )
    subtreelist << subtree + Optional(comma + subtreelist)

    tree = subtree + Word(";").suppress()

    return tree.parseString
#    return hyblabel.parseString

eNewickParser = makeeNewickParser()


if __name__ == '__main__':
    text='((1,,(2)h#LTG1:1.5)x,(h1,3:5.5)y:5.76)r;'
#    text='xx#3'
    expr=eNewickParser(text)
    print(expr)
