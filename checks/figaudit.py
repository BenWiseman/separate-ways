"""Measure a figure's layout instead of judging it.

Two rounds of "the labels look fine to me" shipped a figure whose text ran over the
diagram's own lines, and Ben caught it by reading the PDF. Judgement is not a check. This
records every line a figure draws and every text it places, then reports text sitting on a
line, text sitting on text, and text off the canvas.

Use it like this:

    from figaudit import Panel
    p = Panel(fig, ax)
    p.line(xs, ys, "horizon", color="k")
    p.text(0, 1, "label", ha="center")
    faults = p.audit()          # empty list means clean

`plain_line` and `plain_text` draw without registering, for marks that are meant to sit on
top of something: a number inside a dot on a path, say.

`Panel.selftest()` plants a label in the middle of a registered line and confirms the
audit catches it, because a checker that has never failed is not known to work.
"""
import numpy as np


class Panel:
    def __init__(self, fig, ax, pad=(1.03, 1.10)):
        self.fig, self.ax, self.pad = fig, ax, pad
        self.lines, self.texts = [], []

    # ------------------------------------------------------------------ drawing
    @staticmethod
    def _densify(xs, ys, per_unit=900):
        """Resample a polyline so the audit sees the segments, not just the corners.

        The hexagonal boundary of a Penrose diagram is seven points. Testing a text box
        against those seven caught nothing between them, so a caption resting on a long
        straight edge was reported clean and went to Ben in a PDF. Any check that samples
        a line has to sample all of it.
        """
        xs, ys = np.asarray(xs, float), np.asarray(ys, float)
        if xs.size < 2:
            return xs, ys
        seg = np.hypot(np.diff(xs), np.diff(ys))
        total = float(seg.sum())
        if total <= 0:
            return xs, ys
        n = max(xs.size, int(per_unit * total))
        t = np.concatenate([[0.0], np.cumsum(seg)])
        u = np.linspace(0.0, total, n)
        return np.interp(u, t, xs), np.interp(u, t, ys)

    def line(self, xs, ys, name, **kw):
        self.lines.append((name,) + self._densify(xs, ys))
        return self.ax.plot(xs, ys, **kw)

    def plain_line(self, xs, ys, **kw):
        return self.ax.plot(xs, ys, **kw)

    def text(self, *a, **kw):
        t = self.ax.text(*a, **kw)
        self.texts.append(t)
        return t

    def plain_text(self, *a, **kw):
        return self.ax.text(*a, **kw)

    def annotate(self, *a, **kw):
        t = self.ax.annotate(*a, **kw)
        self.texts.append(t)
        return t

    # ------------------------------------------------------------------- audit
    def _boxes(self, extra=None):
        self.fig.canvas.draw()
        r = self.fig.canvas.get_renderer()
        items = list(self.texts) + ([extra] if extra is not None else [])
        return [(t, t.get_window_extent(renderer=r).expanded(*self.pad)) for t in items]

    def _rect(self, t, r, pad=(1.04, 1.14)):
        """A text's true footprint: centre, width, height, angle, in display coords.

        matplotlib returns the axis-aligned box, which for text at forty-five degrees is
        about (w+h)/sqrt(2) on a side and rejects placements the glyphs clear easily.
        Rotation is about the anchor, and the anchor is the centre only when both
        alignments are "center", so the offset is rotated too.
        """
        import math
        th = math.radians(t.get_rotation())
        keep = t.get_rotation()
        t.set_rotation(0)
        bb = t.get_window_extent(renderer=r)
        t.set_rotation(keep)
        ax_, ay = t.get_transform().transform(t.get_position())
        dx, dy = bb.x0 + bb.width / 2.0 - ax_, bb.y0 + bb.height / 2.0 - ay
        c, sn = math.cos(th), math.sin(th)
        return (ax_ + dx * c - dy * sn, ay + dx * sn + dy * c), \
               bb.width * pad[0], bb.height * pad[1], th

    def audit(self, extra=None):
        import math
        boxes = self._boxes(extra)
        self.fig.canvas.draw()
        r = self.fig.canvas.get_renderer()
        bad = []
        for t, bb in boxes:
            name_t = '"%s"' % " ".join(t.get_text().split())[:30]
            cen, w, h, th = self._rect(t, r)
            for name, xs, ys in self.lines:
                p = self.ax.transData.transform(np.column_stack([xs, ys]))
                du = p - np.asarray(cen)
                u = du[:, 0] * math.cos(th) + du[:, 1] * math.sin(th)
                v = -du[:, 0] * math.sin(th) + du[:, 1] * math.cos(th)
                if np.any((np.abs(u) <= w / 2.0) & (np.abs(v) <= h / 2.0)):
                    bad.append("text %-32s sits on %s" % (name_t, name))
            fb = self.fig.bbox
            if not (bb.x0 >= fb.x0 and bb.x1 <= fb.x1
                    and bb.y0 >= fb.y0 and bb.y1 <= fb.y1):
                bad.append("text %-32s runs off the canvas" % name_t)
        for i in range(len(boxes)):
            for j in range(i + 1, len(boxes)):
                if boxes[i][1].overlaps(boxes[j][1]):
                    bad.append('text "%s" overlaps "%s"'
                               % (" ".join(boxes[i][0].get_text().split())[:26],
                                  " ".join(boxes[j][0].get_text().split())[:26]))
        return sorted(set(bad))

    def place(self, candidates, *a, **kw):
        """Put a text at the first candidate (x, y) that collides with nothing.

        Hand-computed free zones kept being optimistic: a rendered box is wider than the
        arithmetic says, and three rounds of nudging a label by eye is three rounds wasted.
        Searching is cheap and it either finds a clean slot or says there is none.
        """
        t = self.ax.text(candidates[0][0], candidates[0][1], *a, **kw)
        self.texts.append(t)
        for (x, y) in candidates:
            t.set_position((x, y))
            hits = [h for h in self.audit() if self._names(t) in h]
            if not hits:
                return t
        return t                       # leave it at the last try; audit() will report it

    @staticmethod
    def _names(t):
        return '"%s"' % " ".join(t.get_text().split())[:30]

    def selftest(self):
        """Plant a label on a registered line; the audit has to see it."""
        if not self.lines:
            return False
        _, xs, ys = self.lines[0]
        k = len(xs) // 2
        probe = self.ax.text(xs[k], ys[k], "PLANTED", ha="center", va="center",
                             fontsize=9, color="none")
        caught = any("PLANTED" in h for h in self.audit(extra=probe))
        probe.remove()
        return caught
