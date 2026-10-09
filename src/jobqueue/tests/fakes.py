"""Test doubles: a fake ESO archive and a minimal FITS writer."""
import copy
import os

PROGS = {'Io': '60.A-9009(A)', 'Europa': '60.A-9009(B)', 'Ganymede': '60.A-9009(C)', 'Callisto': '60.A-9009(D)'}


class FakeEso:
    """Same methods as eso.EsoArchive. summaries[night][prog] = {'rows', 'by_type', 'last_mod'}."""

    def __init__(self):
        self.summaries = {}
        self.objects = {}
        self.rows = {}           # prog -> list of frame rows (with 'night')
        self.fail = False
        self.calls = []

    def set(self, night, tel, rows, last_mod, by_type=None):
        self.summaries.setdefault(night, {})[PROGS[tel]] = {
            'rows': rows, 'by_type': by_type or {'OBJECT': rows}, 'last_mod': last_mod}

    def night_summary(self, night):
        self.calls.append(('summary', night))
        if self.fail:
            raise RuntimeError('HTTPSConnectionPool: Max retries exceeded')
        return copy.deepcopy(self.summaries.get(night, {}))

    def night_objects(self, prog, night):
        self.calls.append(('objects', prog, night))
        return list(self.objects.get((prog, night), []))

    def frames(self, prog, first, last):
        self.calls.append(('frames', prog, first, last))
        if self.fail:
            raise RuntimeError('TAP HTTP 503')
        return [dict(r) for r in self.rows.get(prog, []) if first <= r['night'] <= last]


def card(key, value):
    if isinstance(value, bool):
        v = '{:>20}'.format('T' if value else 'F')
    elif isinstance(value, (int, float)):
        v = '{:>20}'.format(value)
    else:
        v = "'{}'".format(str(value).replace("'", "''")).ljust(20)
    return '{:8}= {}'.format(key, v).ljust(80)[:80]


def write_fits(path, **keywords):
    cards = [card('SIMPLE', True), card('BITPIX', 16), card('NAXIS', 0)]
    cards += [card(k.replace('_', '-'), v) for k, v in keywords.items()]
    cards.append('END'.ljust(80))
    text = ''.join(cards)
    text += ' ' * (-len(text) % 2880)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'wb') as f:
        f.write(text.encode('ascii'))
    return path
