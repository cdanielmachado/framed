import re
from framed.model.cbmodel import CBModel, CBReaction, Gene
from framed.model.model import Metabolite
from os.path import basename
from itertools import chain
from collections import OrderedDict
from framed.io.sbml import parse_gpr_rule
import warnings

class BiooptParseWarning(Warning):
    pass

def read_cbmodel_from_file(filename, gpr_filename=None):
    parser = BiooptParser()
    model = parser.parse_file(filename)

    if gpr_filename:
        with open(gpr_filename, "r") as gpr_file:
            for line in gpr_file:
                r_id, gpr_string = [col.strip() for col in re.split("\t", line)]
                if r_id in model.reactions and gpr_string:
                    gpr = parse_gpr_rule(gpr_string)
                    for gene in gpr.get_genes():
                        if gene not in model.genes:
                            gene = Gene(elem_id=gene)
                            model.add_gene(gene)

                    model.set_gpr_association(r_id, gpr)
                else:
                    warnings.warn(UserWarning("Reaction {} from GPR was not found in model".format(r_id)))

    return model

class BiooptParser(object):
    def __init__(self, inf=1000):
        # TODO: replace number with float() for performance reasons
        re_number_str = r"(?:[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)|(?:[-+]?(?:[0-9]*\.[0-9]+|[0-9]+))"
        re_bounds_str = "\[\s*(" + re_number_str + r")\s*,\s*(" + re_number_str + r")\s*\]"
        self.re_number = re.compile(re_number_str)
        self.re_constraint = r"(.*)\s*" + re_bounds_str
        self.re_member = re.compile(r"(\(?(" + re_number_str + r") *\)? +)?(.*)")
        self.inf = inf

    def parse_file(self, path):
        """
        :rtype: CBModel
        """
        f = open(path, "r")
        return self.__parse(f.read(), filename=path)

    def parse_reactions_section(self, section_text):
        """
        :rtype: list of :class:`CBReaction`
        """
        return list(r for r, i in self.__parse_reactions_section(section_text))

    def __parse_reactions_section(self, section_text, filename=None, section_start=0, strip_comments=True):
        if not isinstance(section_text, str):
            raise TypeError("Reactions section text is not of type string")

        reactions = self.__parse_section(section_text, lambda x: self.parse_reaction(x, strip_comments=strip_comments), filename=filename, section_start=section_start, strip_comments=strip_comments)

        return reactions

    def parse_reaction_member(self, member_str):
        """
        :rtype: tuple
        """
        if not isinstance(member_str, str):
            raise TypeError("Reaction member string is not of type string")

        member_str = member_str.strip()
        if not len(member_str):
            raise ValueError("Reaction member string is empty")

        m = self.re_member.match(member_str)
        if not m:
            raise SyntaxError("Could not parse reaction member: {0}".format(member_str))

        tmp, coef, name = m.groups()
        coef = float(coef) if coef else 1

        return name, coef

    def parse_reaction_member_list(self, list_str):
        """
        :rtype: list
        """
        if not isinstance(list_str, str):
            raise TypeError("Reaction member list string is not of type string")

        if not len(list_str):
            raise ValueError("Reaction member list string is empty")

        parts = re.split(r"\s+\+\s+", list_str)
        members = (self.parse_reaction_member(s) for s in parts)

        return members

    def parse_reaction(self, line, strip_comments=True):
        """
        :rtype: CBReaction
        """
        if not isinstance(line, str):
            raise TypeError("Reaction line was not a string")

        if strip_comments:
            line, multiline_comment = self.strip_comments(line, False)
        line = line.strip()

        if not len(line):
            raise ValueError("Reaction string is empty")

        sep = ":"
        parts = line.split(sep)

        if len(parts) > 2:
            raise SyntaxError("{0} separator split reaction line into more than two parts [{1}]".format(sep, line))
        if len(parts) < 2:
            raise SyntaxError("Could not split reaction line using {0} separator [{1}]".format(sep, line))

        reaction_name = parts[0].strip()

        d = re.search("(\s+<\->|<\-|\->\s+)", line)
        direction = d.groups()[0]
        parts = parts[1].split(direction)
        direction = direction.strip() #removing white space after splitting into parts

        if len(parts) != 2:
            raise SyntaxError("Reaction doesn't consist of exactly two parts (reactants & products)")

        if direction not in {"<-", "->", "<->"}:
            raise Exception("Unknown direction ({0})".format(direction))

        reactants = ((m_id, -coef) for m_id, coef in self.parse_reaction_member_list(parts[0]))
        products = self.parse_reaction_member_list(parts[1])

        stoichiometry = OrderedDict(chain.from_iterable([reactants, products]))
        reaction = CBReaction(elem_id=reaction_name, stoichiometry=stoichiometry, reversible=direction == "<->")

        return reaction

    # TODO: Add "%" comments functionality
    def strip_comments(self, line, multiline_comment=False):
        short_comment = False
        output_line = ""
        for l in line:
            if multiline_comment and l == "%":
                multiline_comment = False
                continue

            if multiline_comment or short_comment:
                continue

            if l == "%":
                multiline_comment = True
                continue

            if l == "#":
                short_comment = True
                continue

            output_line += l

        return output_line, multiline_comment

    def parse_constraint(self, constraint_text, strip_comments=True):
        """
        :rtype: Bounds
        """
        if strip_comments:
            constraint_text, multiline_comment = self.strip_comments(constraint_text, False)
        constraint_text = constraint_text.strip()

        m = re.match(self.re_constraint, constraint_text)

        if m is None:
            raise SyntaxError("Could parse reaction constraint: {0}".format(constraint_text))

        reaction_name, lb, ub = m.groups()
        reaction_name = reaction_name.strip()

        lb = float(lb)
        lb = None if lb <= -self.inf else lb
        ub = float(ub)
        ub = None if ub >= self.inf else ub

        return reaction_name, lb, ub

    def parse_constraints_section(self, section_text):
        """
        :rtype: list of :class:`Bounds`
        """
        return list(c for c, i in self.__parse_constraints_section(section_text))

    def __parse_constraints_section(self, section_text, filename=None, section_start=0, strip_comments=True):
        if not isinstance(section_text, str):
            raise TypeError("External metabolites section text is not of type string")

        return self.__parse_section(section_text, self.parse_constraint, filename=filename, section_start=section_start, strip_comments=strip_comments)

    def parse_objective_section(self, section_text, filename=None, section_start=0, reactions=None, section_name="-DESIGN OBJECTIVE/-OBJECTIVE", reactions_section_name="-REACTIONS", strip_comments=True):
        if not isinstance(section_text, str):
            raise TypeError("Objective section text is not of type string")

        section_text = re.sub("(\n\r)", "", section_text.strip())
        parts = re.split(r"\s+", section_text)
        if not parts:
            raise SyntaxError("Could not parse objective section line: {0}".format(section_text))

        return parts[0], float(parts[1])


    def parse_external_metabolites_section(self, section_text):
        """
        :rtype: list of :class:`Metabolite`
        """
        return list(e for e, i in self.__parse_external_metabolites_section(section_text))

    def __parse_external_metabolites_section(self, section_text, filename=None, section_start=0, strip_comments=True):
        if not isinstance(section_text, str):
            raise TypeError("External metabolites section text is not of type string")

        return self.__parse_section(section_text, Metabolite, filename=filename, section_start=section_start, strip_comments=strip_comments)

    def find_sections(self, text):
        s = re.compile(r"^-[\w ]+$", re.MULTILINE)
        sections1 = [(m.group(), m.start(), m.end()) for m in re.finditer(s, text)]

        sections2 = dict()
        for i, s in enumerate(sections1):
            name = s[0]
            start = s[2]+1
            end = sections1[i+1][1]-1 if i+1 < len(sections1) else len(text)
            line = text[0:start].count("\n")

            sections2[name] = (start, end, line)

        return sections2

    def __parse_section(self, section_text, method, filename=None, section_start=0, strip_comments=True):
        nl = re.compile("\n\r|\r\n|\n")
        lines = nl.split(section_text)
        comment = False
        for i, line in enumerate(lines):
            if strip_comments:
                line, multiline_comment = self.strip_comments(line, comment)
            line = line.strip()

            # TODO: Return None to have errors with line information anntached
            if not len(line):
                continue

            saved_warnings = []
            with warnings.catch_warnings(record=True) as ws:
                warnings.simplefilter("always")
                result = method(line)

                for w in ws:
                    w_outer = warnings.WarningMessage(message=w.message, category=BiooptParseWarning, filename=filename, lineno=section_start+i+1, line=line)
                    saved_warnings.append(w_outer)

            for w in saved_warnings:
                warnings.warn_explicit(message=w.message, category=w.category, filename=w.filename, lineno=w.lineno)

            yield result, i

    def __find_section(self, text, sections, fun):
        for name in sections.keys():
            if fun(name):
                start, end, line = sections[name]
                return name, text[start:end], line

        return None, None, None

    def parse(self, text):
        """
        :rtype: CBModel
        """
        return self.__parse(text)

    def __parse(self, text, filename=None):
        text = text.replace("\r\n", "\n")
        text = text.replace("\r", "\n")

        re_long = re.compile("%.*?%", re.DOTALL | re.MULTILINE)
        text = re_long.sub("\n", text)

        re_short = re.compile("#.*?\n")
        text = re_short.sub("\n", text)

        sections = self.find_sections(text)

        model = CBModel(basename(filename))
        react_name, react_text, react_line = self.__find_section(text, sections, lambda x: re.search(r"reac", x, re.I))
        const_name, const_text, const_line = self.__find_section(text, sections, lambda x: re.search(r"cons", x, re.I))
        ext_m_name, ext_m_text, ext_m_line = self.__find_section(text, sections, lambda x: re.search(r"ext", x, re.I))
        obj_name, obj_text, obj_line       = self.__find_section(text, sections, lambda x: re.search(r"obj", x, re.I) and not re.search("des", x, re.I))

        if react_text:
            reactions = self.__parse_reactions_section(react_text, filename=filename, section_start=react_line, strip_comments=False)

            for reaction, r_i in reactions:
                for m_id in reaction.stoichiometry:
                    if m_id not in model.metabolites:
                        model.add_metabolite(Metabolite(elem_id=m_id, compartment=None), clear_tmp=True)

                model.add_reaction(reaction, clear_tmp=True)
        else:
            warnings.warn("Could not find '-REACTIONS' section", BiooptParseWarning)

        if const_text:
            for (r_id, lb, ub), i in self.__parse_constraints_section(const_text, filename=filename, section_start=const_line, strip_comments=False):
                if r_id in model.reactions and lb < 0 and not model.reactions[r_id].reversible:
                    warnings.warn_explicit(
                        "Reaction '{0}' from '{1}' has effective bounds not compatible with reaction direction in '{2}' section ({3} : [{4}, {5}])".format(r_id, const_name, react_name, "<->" if reactions[r_id].reversibl else "->", lb, ub),
                        BiooptParseWarning, filename=filename, lineno=const_line+i+1)

                if r_id in model.reactions:
                    model.reactions[r_id].lb = lb
                    model.reactions[r_id].ub = ub
                elif react_text:
                    warnings.warn_explicit(
                        "Reaction '{0}' from '{1}' section is not present in '{2}' section".format(r_id, const_name, react_name),
                        BiooptParseWarning, filename=filename, lineno=const_line+i+1)
        else:
            warnings.warn("Could not find '-CONSTRAINS' section", BiooptParseWarning)

        if ext_m_text:
            met2rxn = model.metabolite_reaction_lookup()
            for m_external, i in self.__parse_external_metabolites_section(ext_m_text, filename=filename, section_start=ext_m_line, strip_comments=False):
                if m_external.id in model.metabolites:
                    model.metabolites[m_external.id].boundary = True
                    for external_r_id, coef in met2rxn[m_external.id].iteritems():
                        model.reactions[external_r_id].is_exchange = True
                        del model.reactions[external_r_id].stoichiometry[m_external.id]
                elif react_text:
                    warnings.warn_explicit(
                        "Metabolite '{0}' from '{1}' section is not present in any reaction from '{2}' section".format(m.id, ext_m_name, react_name),
                        BiooptParseWarning, filename=filename, lineno=ext_m_line+i+1)
        else:
            warnings.warn("Could not find '-EXTERNAL METABOLITES' section", BiooptParseWarning)

        if obj_text:
            r_id, coef = self.parse_objective_section(obj_text, section_name=obj_name, reactions_section_name=react_name, filename=filename, section_start=obj_line, reactions=reactions, strip_comments=False)
            model.reactions[r_id].objective = coef
        else:
            warnings.warn("Could not find '-OBJECTIVE' section", BiooptParseWarning)


        return model