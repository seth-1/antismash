# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Annotations for secondary metabolites """

import re
from typing import Any, Iterable, Iterator, List, Sequence, Set, Union


def _parse_format(fmt: str, data: str) -> Sequence[str]:
    """ Reverse of str.format(), pulls values from an input string that match
        positions of {} in a format string. Raises a ValueError if the match
        cannot be found.
    """
    safe = fmt.replace('(', r'\(').replace(')', r'\)')
    regex = "^{}$".format(safe.replace("{}", "(.+?)"))
    res = re.search(regex, data)
    if res is None:
        raise ValueError("Could not match format %r to input %r" % (fmt, data))
    return res.groups()


class SecMetQualifier(list):
    """ A qualifier for tracking various secondary metabolite information about
        a CDS.

        Can be directly used as a qualifier for BioPython's SeqFeature.
    """
    class Domain:
        """ A simple container for the information needed to create a domain """
        qualifier_label = "{} (E-value: {}, bitscore: {}, seeds: {}, tool: {})"

        def __init__(self, name: str, evalue: float, bitscore: float, nseeds: int,
                     tool: str) -> None:
            self.query_id = str(name)
            self.evalue = float(evalue)
            self.bitscore = float(bitscore)
            self.nseeds = int(nseeds)
            self.tool = str(tool)

        def __repr__(self) -> str:
            return str(self)

        def __str__(self) -> str:
            return self.qualifier_label.format(self.query_id, self.evalue,
                                               self.bitscore, self.nseeds, self.tool)

        def to_json(self) -> List[Union[str, float, int]]:
            """ Constructs a JSON-friendly representation of a Domain """
            return [self.query_id, self.evalue, self.bitscore, self.nseeds, self.tool]

        @classmethod
        def from_string(cls, line: str) -> "SecMetQualifier.Domain":
            """ Rebuilds a Domain from a string (e.g. from a genbank file) """
            return cls.from_json(_parse_format(cls.qualifier_label, line))

        @classmethod
        def from_json(cls, json: Sequence[Union[str, float]]) -> "SecMetQualifier.Domain":
            """ Rebuilds a Domain from a JSON representation """
            assert len(json) == 5, json
            return cls(str(json[0]), float(json[1]), float(json[2]), int(json[3]), str(json[4]))

    def __init__(self, products: Set[str], domains: List["SecMetQualifier.Domain"]) -> None:
        self._domains = domains
        self.domain_ids = []  # type: List[str]
        self.unique_domain_ids = set()  # type: Set[str]
        for domain in self._domains:
            assert isinstance(domain, SecMetQualifier.Domain)
            assert domain.query_id not in self.unique_domain_ids, "domains were duplicated"
            self.unique_domain_ids.add(domain.query_id)
            self.domain_ids.append(domain.query_id)
        self._products = set()  # type: Set[str]
        self.add_products(products)
        self.kind = "biosynthetic"
        super().__init__()

    def __iter__(self) -> Iterator[str]:
        yield "Type: %s" % self.clustertype
        yield "; ".join(map(str, self._domains))
        yield "Kind: %s" % self.kind

    def append(self, _item: Any) -> None:
        raise NotImplementedError("Appending to this list won't work")

    def extend(self, _items: Iterable[Any]) -> None:
        raise NotImplementedError("Extending this list won't work")

    def add_products(self, products: Set[str]) -> None:
        """ Adds one or more products to the qualifier """
        assert isinstance(products, set), type(products)
        for product in products:
            assert isinstance(product, str) and "-" not in product, product
        self._products.update(products)

    def add_domains(self, domains: List["SecMetQualifier.Domain"]) -> None:
        """ Add a group of Domains to the the qualifier """
        unique = []
        for domain in domains:
            assert isinstance(domain, SecMetQualifier.Domain)
            if domain.query_id in self.unique_domain_ids:
                continue  # no sense keeping duplicates
            self.unique_domain_ids.add(domain.query_id)
            unique.append(domain)
        self._domains.extend(unique)

    @property
    def domains(self) -> List["SecMetQualifier.Domain"]:
        """ A list of domains stored in the qualifier"""
        return list(self._domains)

    @property
    def products(self) -> List[str]:
        """ A list of all products a feature is involved in"""
        return sorted(self._products)

    @property
    def clustertype(self) -> str:
        """ A string hypen-separated products """
        return "-".join(sorted(self.products))

    @staticmethod
    def from_biopython(qualifier: List[str]) -> "SecMetQualifier":
        """ Converts a BioPython style qualifier into a SecMetQualifier. """
        domains = []
        products = set()  # type: Set[str]
        kind = "biosynthetic"
        if len(qualifier) != 3:
            raise ValueError("Cannot parse qualifier: %s" % qualifier)
        for value in qualifier:
            if value.startswith("Type: "):
                products = set(value.split("Type: ", 1)[1].split("-"))
            elif value.startswith("Kind: "):
                kind = value.split("Kind: ", 1)[1]
                assert kind == "biosynthetic", kind  # since it's the only kind we have
            else:
                domain_strings = value.split("; ")
                for domain_string in domain_strings:
                    domains.append(SecMetQualifier.Domain.from_string(domain_string))
        if not (domains and products and kind):
            raise ValueError("Cannot parse qualifier: %s" % qualifier)
        return SecMetQualifier(products, domains)

    def __len__(self) -> int:
        return 3
