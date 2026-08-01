// Self-contained HTML run report.
//
// The page is meant to be read by someone who has never seen the source: every
// stage carries a sentence saying what it does before the numbers that say what
// it did. Charts are emitted as inline SVG built here in C++ -- no scripts are
// fetched, no fonts are loaded, nothing on the page touches the network.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "report.h"

namespace ts {
namespace {

// ---------------------------------------------------------------- formatting

std::string htmlEscape(const std::string& in) {
    std::string out;
    out.reserve(in.size() + 8);
    for (char c : in) {
        switch (c) {
            case '&':  out += "&amp;";  break;
            case '<':  out += "&lt;";   break;
            case '>':  out += "&gt;";   break;
            case '"':  out += "&quot;"; break;
            case '\'': out += "&#39;";  break;
            default:   out += c;
        }
    }
    return out;
}

std::string groupDigits(const std::string& digits) {
    std::string out;
    int run = 0;
    for (size_t i = digits.size(); i-- > 0;) {
        out += digits[i];
        if (++run % 3 == 0 && i > 0) out += ',';
    }
    std::reverse(out.begin(), out.end());
    return out;
}

std::string fmtNum(double v, int decimals) {
    if (!std::isfinite(v)) return "n/a";
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", decimals, v);
    std::string s(buf);
    const bool neg = !s.empty() && s[0] == '-';
    if (neg) s.erase(0, 1);
    std::string ip = s, fp;
    const size_t dot = s.find('.');
    if (dot != std::string::npos) { ip = s.substr(0, dot); fp = s.substr(dot); }
    std::string r = groupDigits(ip) + fp;
    if (neg && r != "0") r.insert(r.begin(), '-');
    return r;
}

std::string fmtInt(double v) { return fmtNum(v, 0); }

// Every ratio in this file goes through here, so a zero denominator can never
// reach the page as NaN or inf.
double ratio(double num, double den) { return den > 0 ? num / den : 0.0; }

// Same guard, but for spans that are legitimately negative (SVG y grows down).
double frac(double num, double den) { return den != 0 ? num / den : 0.0; }

std::string fmtPct(double num, double den, int decimals = 1) {
    if (!(den > 0)) return "n/a";
    return fmtNum(100.0 * num / den, decimals) + "%";
}

std::string fmtPctValue(double percent, int decimals = 1) {
    return fmtNum(percent, decimals) + "%";
}

std::string fmtBp(double bp) {
    if (bp >= 1e9) return fmtNum(bp / 1e9, 2) + " Gb";
    if (bp >= 1e6) return fmtNum(bp / 1e6, 2) + " Mb";
    if (bp >= 1e3) return fmtNum(bp / 1e3, 1) + " kb";
    return fmtInt(bp) + " bp";
}

std::string fmtBytes(double b) {
    if (b >= 1024.0 * 1024 * 1024) return fmtNum(b / (1024.0 * 1024 * 1024), 2) + " GiB";
    if (b >= 1024.0 * 1024)        return fmtNum(b / (1024.0 * 1024), 1) + " MiB";
    if (b >= 1024.0)               return fmtNum(b / 1024.0, 1) + " KiB";
    return fmtInt(b) + " B";
}

std::string fmtDuration(double sec) {
    if (!std::isfinite(sec) || sec < 0) return "n/a";
    if (sec < 1.0)   return fmtNum(sec * 1000.0, 0) + " ms";
    if (sec < 90.0)  return fmtNum(sec, 2) + " s";
    const long total = static_cast<long>(sec + 0.5);
    const long h = total / 3600, m = (total % 3600) / 60, s = total % 60;
    char buf[64];
    if (h > 0) std::snprintf(buf, sizeof(buf), "%ld h %02ld m %02ld s", h, m, s);
    else       std::snprintf(buf, sizeof(buf), "%ld m %02ld s", m, s);
    return buf;
}

// Compact axis-tick text: 12,300 -> "12k".
std::string fmtTick(double v) {
    const double a = std::fabs(v);
    if (a >= 1e9) return fmtNum(v / 1e9, a < 1e10 ? 1 : 0) + "G";
    if (a >= 1e6) return fmtNum(v / 1e6, a < 1e7 ? 1 : 0) + "M";
    if (a >= 1e3) return fmtNum(v / 1e3, a < 1e4 ? 1 : 0) + "k";
    if (a >= 10)  return fmtNum(v, 0);
    if (a >= 1)   return fmtNum(v, v == std::floor(v) ? 0 : 1);
    if (a == 0)   return "0";
    return fmtNum(v, 2);
}

// Short numeric text for SVG geometry attributes.
std::string n2(double v) {
    if (!std::isfinite(v)) v = 0;
    char buf[48];
    std::snprintf(buf, sizeof(buf), "%.2f", v);
    std::string s(buf);
    if (s.find('.') != std::string::npos) {
        while (!s.empty() && s.back() == '0') s.pop_back();
        if (!s.empty() && s.back() == '.') s.pop_back();
    }
    if (s == "-0" || s.empty()) s = "0";
    return s;
}

std::string at(const char* name, double v) {
    return std::string(" ") + name + "=\"" + n2(v) + "\"";
}
std::string at(const char* name, const std::string& v) {
    return std::string(" ") + name + "=\"" + v + "\"";
}

double clamp01(double v) { return v < 0 ? 0 : (v > 1 ? 1 : v); }

// Position of v on a log scale running from lo (0) to hi (1).
double logSpan(double v, double lo, double hi) {
    if (!(v > 0) || !(lo > 0) || !(hi > lo)) return 0;
    return clamp01((std::log10(v) - std::log10(lo)) / (std::log10(hi) - std::log10(lo)));
}

// ---------------------------------------------------------------- chart core

struct Axis {
    double lo = 0, hi = 1;
    double p0 = 0, p1 = 1;   // pixel of lo, pixel of hi
    bool log = false;

    double t(double v) const { return log ? std::log10(v > 1e-9 ? v : 1e-9) : v; }
    double px(double v) const {
        const double a = t(lo), b = t(hi);
        const double f = (b - a) > 0 ? (t(v) - a) / (b - a) : 0.0;
        return p0 + f * (p1 - p0);
    }
    // Pads both ends, for an axis that tracks a series rather than starting at
    // zero: without it a line hugging its own min or max sits on the frame.
    void fitRange(double dataLo, double dataHi, double padFrac) {
        log = false;
        if (!std::isfinite(dataLo) || !std::isfinite(dataHi)) { lo = 0; hi = 1; return; }
        if (!(dataHi > dataLo)) {
            const double span = std::fabs(dataLo) > 0 ? std::fabs(dataLo) * 0.1 : 1.0;
            lo = dataLo - span;
            hi = dataLo + span;
            return;
        }
        const double pad = (dataHi - dataLo) * padFrac;
        lo = dataLo - pad;
        hi = dataHi + pad;
    }

    void fit(double dataLo, double dataHi, bool isLog, double padFrac = 0.05) {
        log = isLog;
        if (!std::isfinite(dataLo) || !std::isfinite(dataHi)) { lo = 0; hi = 1; return; }
        if (isLog) {
            lo = dataLo > 0 ? dataLo : 1.0;
            hi = dataHi > lo ? dataHi : lo * 10.0;
            lo = std::pow(10.0, std::floor(std::log10(lo)));
            hi = std::pow(10.0, std::ceil(std::log10(hi)));
        } else {
            lo = dataLo;
            hi = dataHi;
            if (!(hi > lo)) hi = lo + (std::fabs(lo) > 0 ? std::fabs(lo) : 1.0);
            const double pad = (hi - lo) * padFrac;
            hi += pad;
        }
    }
};

std::vector<double> linearTicks(double lo, double hi, int target) {
    std::vector<double> out;
    if (!(hi > lo) || target < 1) { out.push_back(lo); return out; }
    const double raw = (hi - lo) / target;
    const double mag = std::pow(10.0, std::floor(std::log10(raw)));
    const double norm = raw / mag;
    double step = mag;
    if (norm > 5)      step = 10 * mag;
    else if (norm > 2) step = 5 * mag;
    else if (norm > 1) step = 2 * mag;
    if (!(step > 0)) { out.push_back(lo); return out; }
    const double start = std::ceil(lo / step) * step;
    for (double v = start; v <= hi + step * 1e-6 && out.size() < 40; v += step) {
        out.push_back(std::fabs(v) < step * 1e-9 ? 0.0 : v);
    }
    if (out.empty()) out.push_back(lo);
    return out;
}

std::vector<double> logTicks(double lo, double hi) {
    std::vector<double> out;
    if (!(lo > 0) || !(hi > lo)) { out.push_back(lo > 0 ? lo : 1); return out; }
    const int a = static_cast<int>(std::floor(std::log10(lo)));
    const int b = static_cast<int>(std::ceil(std::log10(hi)));
    for (int e = a; e <= b && out.size() < 20; ++e) {
        const double v = std::pow(10.0, e);
        if (v >= lo * 0.999 && v <= hi * 1.001) out.push_back(v);
    }
    if (out.empty()) out.push_back(lo);
    return out;
}

struct ChartOpts {
    std::string title;       // already escaped
    std::string subtitle;    // already escaped
    std::string xlab;
    std::string ylab;
    std::string y2lab;
    std::string y3lab;
    double w = 820, h = 330;
    double mL = 64, mR = 18, mT = 12, mB = 48;
    double minWidth = 500;   // below this the chart's own container scrolls
};

// The SVG is stretched to its container width, so the viewBox has to be close
// to the width the chart will actually get or every label scales with it. This
// is the geometry for a chart sharing a row with another one.
ChartOpts narrowChart() {
    ChartOpts o;
    o.w = 520;
    o.h = 300;
    o.mL = 56;
    o.mR = 14;
    o.mT = 12;
    o.mB = 46;
    o.minWidth = 360;
    return o;
}

// Draws the frame (grid, axes, ticks, labels) around marks supplied by the
// caller. Up to two right-hand tick columns exist so a chart can overlay two
// series that share no unit with the left axis.
struct Frame {
    explicit Frame(ChartOpts opts) : o(std::move(opts)) {
        x.p0 = o.mL;        x.p1 = o.w - o.mR;
        y.p0 = o.h - o.mB;  y.p1 = o.mT;
        y2 = y3 = y;
    }

    ChartOpts o;
    Axis x, y, y2, y3;

    std::string underlay;   // shading, drawn beneath the grid
    std::string marks;      // the data itself
    std::string overlay;    // reference lines and labels, drawn last
    std::string legend;     // HTML, sits under the plot

    std::vector<double> xTicks, yTicks, y2Ticks, y3Ticks;
    std::vector<std::string> xTickLabels, yTickLabels, y2TickLabels, y3TickLabels;
    int y2Color = 5, y3Color = 6;
    bool xGrid = false;
    bool rotateXTicks = false;

    double left()   const { return o.mL; }
    double right()  const { return o.w - o.mR; }
    double top()    const { return o.mT; }
    double bottom() const { return o.h - o.mB; }

    std::string render() const;
};

std::string Frame::render() const {
    const double L = left(), R = right(), T = top(), B = bottom();
    std::string s;
    s += "<figure class=\"chart\">";
    if (!o.title.empty()) s += "<figcaption class=\"chart-title\">" + o.title + "</figcaption>";
    if (!o.subtitle.empty()) s += "<p class=\"chart-sub\">" + o.subtitle + "</p>";
    s += "<div class=\"chart-wrap\"><svg class=\"cv\" role=\"img\"";
    s += at("viewBox", "0 0 " + n2(o.w) + " " + n2(o.h));
    s += at("style", "min-width:" + n2(o.minWidth) + "px");
    s += ">";
    if (!o.title.empty()) s += "<title>" + o.title + "</title>";

    s += underlay;

    for (double v : yTicks) {
        const double py = y.px(v);
        s += "<line class=\"grid\"" + at("x1", L) + at("y1", py) + at("x2", R) + at("y2", py) + "/>";
    }
    if (xGrid) {
        for (double v : xTicks) {
            const double pxv = x.px(v);
            s += "<line class=\"grid\"" + at("x1", pxv) + at("y1", T) + at("x2", pxv) + at("y2", B) + "/>";
        }
    }

    s += "<line class=\"axis\"" + at("x1", L) + at("y1", B) + at("x2", R) + at("y2", B) + "/>";
    s += "<line class=\"axis\"" + at("x1", L) + at("y1", T) + at("x2", L) + at("y2", B) + "/>";

    s += marks;
    s += overlay;

    for (size_t i = 0; i < yTicks.size(); ++i) {
        const double py = y.px(yTicks[i]);
        const std::string lab = i < yTickLabels.size() ? yTickLabels[i] : fmtTick(yTicks[i]);
        s += "<text class=\"tick\"" + at("x", L - 7) + at("y", py + 3.6) +
             " text-anchor=\"end\">" + lab + "</text>";
    }
    for (size_t i = 0; i < xTicks.size(); ++i) {
        const double pxv = x.px(xTicks[i]);
        const std::string lab = i < xTickLabels.size() ? xTickLabels[i] : fmtTick(xTicks[i]);
        if (rotateXTicks) {
            s += "<text class=\"tick\" text-anchor=\"end\"" +
                 at("transform", "translate(" + n2(pxv + 3) + "," + n2(B + 13) + ") rotate(-40)") +
                 ">" + lab + "</text>";
        } else {
            s += "<text class=\"tick\"" + at("x", pxv) + at("y", B + 15) +
                 " text-anchor=\"middle\">" + lab + "</text>";
        }
    }
    for (size_t i = 0; i < y2Ticks.size(); ++i) {
        const double py = y2.px(y2Ticks[i]);
        const std::string lab = i < y2TickLabels.size() ? y2TickLabels[i] : fmtTick(y2Ticks[i]);
        s += "<text class=\"tick t" + std::to_string(y2Color) + "\"" + at("x", R + 6) +
             at("y", py + 3.6) + " text-anchor=\"start\">" + lab + "</text>";
    }
    for (size_t i = 0; i < y3Ticks.size(); ++i) {
        const double py = y3.px(y3Ticks[i]);
        const std::string lab = i < y3TickLabels.size() ? y3TickLabels[i] : fmtTick(y3Ticks[i]);
        s += "<text class=\"tick t" + std::to_string(y3Color) + "\"" + at("x", R + 48) +
             at("y", py + 3.6) + " text-anchor=\"start\">" + lab + "</text>";
    }

    if (!o.xlab.empty()) {
        s += "<text class=\"axlab\"" + at("x", (L + R) / 2) + at("y", o.h - 6) +
             " text-anchor=\"middle\">" + o.xlab + "</text>";
    }
    if (!o.ylab.empty()) {
        s += "<text class=\"axlab\" text-anchor=\"middle\"" +
             at("transform", "translate(12," + n2((T + B) / 2) + ") rotate(-90)") +
             ">" + o.ylab + "</text>";
    }
    if (!o.y2lab.empty()) {
        s += "<text class=\"axlab t" + std::to_string(y2Color) + "\" text-anchor=\"middle\"" +
             at("transform", "translate(" + n2(o.w - 8) + "," + n2((T + B) / 2) + ") rotate(90)") +
             ">" + o.y2lab + "</text>";
    }
    if (!o.y3lab.empty()) {
        s += "<text class=\"axlab t" + std::to_string(y3Color) + "\" text-anchor=\"middle\"" +
             at("transform", "translate(" + n2(o.w - 26) + "," + n2((T + B) / 2) + ") rotate(90)") +
             ">" + o.y3lab + "</text>";
    }

    s += "</svg></div>";
    if (!legend.empty()) s += "<div class=\"legend\">" + legend + "</div>";
    s += "</figure>";
    return s;
}

std::string legendSwatch(int color, const std::string& name, bool line = false) {
    return "<span class=\"lg\"><i class=\"sw" + std::string(line ? " swline" : "") +
           " bg" + std::to_string(color) + "\"></i>" + name + "</span>";
}

std::string noDataChart(const std::string& title, const std::string& why) {
    return "<figure class=\"chart\"><figcaption class=\"chart-title\">" + title +
           "</figcaption><p class=\"nodata\">" + why + "</p></figure>";
}

// ---------------------------------------------------------------- primitives

std::string svgRect(double xx, double yy, double w, double h, const std::string& cls,
                    const std::string& tip) {
    if (w < 0) { xx += w; w = -w; }
    if (h < 0) { yy += h; h = -h; }
    std::string s = "<rect" + at("x", xx) + at("y", yy) + at("width", w) + at("height", h) +
                    at("class", cls);
    if (!tip.empty()) s += at("data-tip", tip);
    s += "/>";
    return s;
}

std::string svgPolyline(const std::vector<double>& px, const std::vector<double>& py,
                        const std::string& cls) {
    if (px.size() < 2) return std::string();
    std::string pts;
    for (size_t i = 0; i < px.size() && i < py.size(); ++i) {
        if (i) pts += ' ';
        pts += n2(px[i]) + "," + n2(py[i]);
    }
    return "<polyline" + at("class", cls) + at("points", pts) + "/>";
}

std::string svgVMarker(const Frame& f, double xv, const std::string& label, int color,
                       double labelDy, const std::string& tip) {
    const double pxv = f.x.px(xv);
    if (pxv < f.left() - 1 || pxv > f.right() + 1) return std::string();
    std::string s = "<line class=\"vmark k" + std::to_string(color) + "\"" + at("x1", pxv) +
                    at("y1", f.top()) + at("x2", pxv) + at("y2", f.bottom());
    if (!tip.empty()) s += at("data-tip", tip);
    s += "/>";
    const bool flip = pxv > (f.left() + f.right()) * 0.62;
    s += "<text class=\"vlab t" + std::to_string(color) + "\"" + at("x", pxv + (flip ? -5 : 5)) +
         at("y", f.top() + labelDy) + at("text-anchor", flip ? "end" : "start") + ">" + label +
         "</text>";
    return s;
}

// ------------------------------------------------------------- chart helpers

struct HistMarker {
    double x = 0;
    std::string label;
    int color = 5;
    double labelDy = 12;
    std::string tip;
};

// Bars for a dense integer-indexed distribution. Above ~260 buckets the bars
// stop being individually visible, so the series is drawn as a filled area with
// coarser hover targets instead.
std::string svgHistogram(ChartOpts opts, const std::vector<double>& values, double x0,
                         double xStep, bool logY, int color,
                         const std::string& xUnit, const std::string& yUnit,
                         const std::vector<HistMarker>& markers = {},
                         double shadeFrom = 0, double shadeTo = 0,
                         const std::string& shadeTip = std::string()) {
    if (values.empty()) return noDataChart(opts.title, "No data was recorded for this chart.");

    double maxV = 0;
    for (double v : values) maxV = std::max(maxV, v);
    if (!(maxV > 0)) return noDataChart(opts.title, "Every bucket in this distribution is empty.");

    Frame f(std::move(opts));
    const size_t n = values.size();
    f.x.fit(x0, x0 + xStep * static_cast<double>(n), false, 0.0);
    f.x.hi = x0 + xStep * static_cast<double>(n);
    if (logY) {
        // Baseline slightly below 1 so single-count buckets still show a bar.
        f.y.log = true;
        f.y.lo = 0.6;
        f.y.hi = std::pow(10.0, std::ceil(std::log10(maxV)));
        if (!(f.y.hi > 1)) f.y.hi = 10;
        f.yTicks = logTicks(1.0, f.y.hi);
    } else {
        f.y.fit(0, maxV, false, 0.08);
        f.yTicks = linearTicks(0, f.y.hi, 5);
    }
    f.xTicks = linearTicks(f.x.lo, f.x.hi, 10);

    if (shadeTo > shadeFrom) {
        const double a = f.x.px(std::max(shadeFrom, f.x.lo));
        const double b = f.x.px(std::min(shadeTo, f.x.hi));
        f.underlay += svgRect(a, f.top(), b - a, f.bottom() - f.top(), "shade", shadeTip);
    }

    const double base = f.y.px(logY ? 0.6 : 0.0);
    const std::string fillCls = "bar f" + std::to_string(color);

    if (n <= 260) {
        const double bw = (f.right() - f.left()) / static_cast<double>(n);
        for (size_t i = 0; i < n; ++i) {
            if (!(values[i] > 0)) continue;
            const double xv = x0 + xStep * static_cast<double>(i);
            const double px0 = f.x.px(xv);
            const double top = f.y.px(values[i]);
            const std::string tip = fmtTick(xv) + " " + xUnit + " \xE2\x80\x94 " +
                                    fmtInt(values[i]) + " " + yUnit;
            f.marks += svgRect(px0, top, std::max(bw - 0.6, 0.8), base - top, fillCls, tip);
        }
    } else {
        std::vector<double> pxs, pys;
        pxs.push_back(f.x.px(x0));
        pys.push_back(base);
        for (size_t i = 0; i < n; ++i) {
            const double xv = x0 + xStep * static_cast<double>(i);
            pxs.push_back(f.x.px(xv));
            pys.push_back(f.y.px(values[i] > 0 ? values[i] : (logY ? 0.6 : 0.0)));
        }
        pxs.push_back(f.x.px(x0 + xStep * static_cast<double>(n - 1)));
        pys.push_back(base);
        std::string pts;
        for (size_t i = 0; i < pxs.size(); ++i) {
            if (i) pts += ' ';
            pts += n2(pxs[i]) + "," + n2(pys[i]);
        }
        f.marks += "<polygon" + at("class", "area a" + std::to_string(color)) + at("points", pts) + "/>";

        // Aggregated hover targets, roughly one every 6 px.
        const double plotW = f.right() - f.left();
        const size_t groups = static_cast<size_t>(std::max(1.0, plotW / 6.0));
        const size_t per = std::max<size_t>(1, n / groups);
        for (size_t i = 0; i < n; i += per) {
            const size_t j = std::min(n, i + per);
            double sum = 0, peak = 0;
            size_t peakAt = i;
            for (size_t q = i; q < j; ++q) {
                sum += values[q];
                if (values[q] > peak) { peak = values[q]; peakAt = q; }
            }
            if (!(sum > 0)) continue;
            const double xa = x0 + xStep * static_cast<double>(i);
            const double xb = x0 + xStep * static_cast<double>(j - 1);
            const std::string tip = fmtTick(xa) + "\xE2\x80\x93" + fmtTick(xb) + " " + xUnit +
                                    " \xE2\x80\x94 " + fmtInt(sum) + " " + yUnit + ", peak at " +
                                    fmtTick(x0 + xStep * static_cast<double>(peakAt)) + " " + xUnit;
            f.marks += svgRect(f.x.px(xa), f.top(), std::max(f.x.px(xb) - f.x.px(xa), 1.0),
                               f.bottom() - f.top(), "hit", tip);
        }
    }

    for (const HistMarker& m : markers) {
        f.overlay += svgVMarker(f, m.x, m.label, m.color, m.labelDy, m.tip);
    }
    return f.render();
}

std::string svgBarChart(ChartOpts opts, const std::vector<std::string>& labels,
                        const std::vector<double>& values, bool logY, int color,
                        const std::vector<std::string>& tips) {
    if (values.empty()) return noDataChart(opts.title, "No data was recorded for this chart.");
    double maxV = 0;
    for (double v : values) maxV = std::max(maxV, v);
    if (!(maxV > 0)) return noDataChart(opts.title, "Every bucket in this distribution is empty.");

    Frame f(std::move(opts));
    f.rotateXTicks = true;
    const size_t n = values.size();
    f.x.lo = 0;
    f.x.hi = static_cast<double>(n);
    if (logY) {
        f.y.log = true;
        f.y.lo = 0.6;
        f.y.hi = std::pow(10.0, std::ceil(std::log10(maxV)));
        if (!(f.y.hi > 1)) f.y.hi = 10;
        f.yTicks = logTicks(1.0, f.y.hi);
    } else {
        f.y.fit(0, maxV, false, 0.08);
        f.yTicks = linearTicks(0, f.y.hi, 5);
    }
    const double base = f.y.px(logY ? 0.6 : 0.0);
    const double bw = (f.right() - f.left()) / static_cast<double>(n);
    for (size_t i = 0; i < n; ++i) {
        f.xTicks.push_back(static_cast<double>(i) + 0.5);
        f.xTickLabels.push_back(i < labels.size() ? labels[i] : std::string());
        if (!(values[i] > 0)) continue;
        const double top = f.y.px(values[i]);
        f.marks += svgRect(f.x.px(static_cast<double>(i)) + bw * 0.12, top, bw * 0.76, base - top,
                           "bar f" + std::to_string(color), i < tips.size() ? tips[i] : std::string());
    }
    return f.render();
}

struct LineSeries {
    std::string name;
    std::vector<double> y;
    std::vector<std::string> tips;
    int color = 1;
    int axis = 0;      // 0 = left, 1 = right column 1, 2 = right column 2
    bool step = false;
    bool points = true;
};

std::string svgLineChart(ChartOpts opts, const std::vector<double>& xs,
                         const std::vector<LineSeries>& series,
                         const std::vector<std::string>& xTickLabels,
                         bool logY = false) {
    if (xs.empty() || series.empty()) {
        return noDataChart(opts.title, "No data was recorded for this chart.");
    }
    Frame f(std::move(opts));

    double xlo = xs[0], xhi = xs[0];
    for (double v : xs) { xlo = std::min(xlo, v); xhi = std::max(xhi, v); }
    f.x.lo = xlo;
    f.x.hi = xhi > xlo ? xhi : xlo + 1;
    const double xpad = (f.x.hi - f.x.lo) * 0.04;
    f.x.lo -= xpad;
    f.x.hi += xpad;

    double lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
    bool has[3] = {false, false, false};
    for (const LineSeries& sr : series) {
        const int ax = sr.axis < 0 || sr.axis > 2 ? 0 : sr.axis;
        for (double v : sr.y) {
            if (!std::isfinite(v)) continue;
            if (!has[ax]) { lo[ax] = hi[ax] = v; has[ax] = true; }
            lo[ax] = std::min(lo[ax], v);
            hi[ax] = std::max(hi[ax], v);
        }
    }
    Axis* axes[3] = {&f.y, &f.y2, &f.y3};
    for (int i = 0; i < 3; ++i) {
        if (!has[i]) continue;
        if (logY && i == 0) {
            axes[i]->fit(lo[i] > 0 ? lo[i] : 1, hi[i], true);
        } else {
            axes[i]->fit(std::min(0.0, lo[i]), hi[i], false, 0.10);
        }
    }
    f.yTicks = f.y.log ? logTicks(f.y.lo, f.y.hi) : linearTicks(f.y.lo, f.y.hi, 5);
    // The right-hand columns reuse the left gridline positions, so the chart
    // keeps a single set of horizontal rules however many units it shows.
    for (int i = 1; i < 3; ++i) {
        if (!has[i]) continue;
        std::vector<double>& dst = (i == 1) ? f.y2Ticks : f.y3Ticks;
        for (double t : f.yTicks) {
            const double g = frac(f.y.px(t) - f.y.p0, f.y.p1 - f.y.p0);
            dst.push_back(axes[i]->lo + g * (axes[i]->hi - axes[i]->lo));
        }
    }
    for (const LineSeries& sr : series) {
        if (sr.axis == 1) f.y2Color = sr.color;
        if (sr.axis == 2) f.y3Color = sr.color;
    }

    f.xTicks = xs;
    f.xTickLabels = xTickLabels;
    if (f.xTickLabels.empty()) for (double v : xs) f.xTickLabels.push_back(fmtTick(v));
    f.xGrid = true;

    for (const LineSeries& sr : series) {
        const int ax = sr.axis < 0 || sr.axis > 2 ? 0 : sr.axis;
        const Axis& ya = *axes[ax];
        std::vector<double> px, py;
        for (size_t i = 0; i < xs.size() && i < sr.y.size(); ++i) {
            if (!std::isfinite(sr.y[i])) continue;
            if (sr.step && !px.empty()) { px.push_back(f.x.px(xs[i])); py.push_back(py.back()); }
            px.push_back(f.x.px(xs[i]));
            py.push_back(ya.px(sr.y[i]));
        }
        f.marks += svgPolyline(px, py, "ln k" + std::to_string(sr.color));
        if (sr.points) {
            for (size_t i = 0; i < xs.size() && i < sr.y.size(); ++i) {
                if (!std::isfinite(sr.y[i])) continue;
                std::string s = "<circle" + at("cx", f.x.px(xs[i])) + at("cy", ya.px(sr.y[i])) +
                                at("r", 3.4) + at("class", "pt f" + std::to_string(sr.color));
                if (i < sr.tips.size() && !sr.tips[i].empty()) s += at("data-tip", sr.tips[i]);
                s += "/>";
                f.marks += s;
            }
        }
        f.legend += legendSwatch(sr.color, sr.name, true);
    }
    return f.render();
}

struct StackSeries {
    std::string name;
    std::vector<double> v;
    int color = 1;
};

std::string svgStackedBar(ChartOpts opts, const std::vector<std::string>& groupLabels,
                          const std::vector<StackSeries>& stack,
                          const std::vector<LineSeries>& overlays,
                          const std::string& stackUnit) {
    if (groupLabels.empty() || stack.empty()) {
        return noDataChart(opts.title, "No data was recorded for this chart.");
    }
    const size_t n = groupLabels.size();
    std::vector<double> totals(n, 0.0);
    for (const StackSeries& sr : stack) {
        for (size_t i = 0; i < n && i < sr.v.size(); ++i) totals[i] += sr.v[i];
    }
    double maxTotal = 0;
    for (double v : totals) maxTotal = std::max(maxTotal, v);

    Frame f(std::move(opts));
    f.x.lo = 0;
    f.x.hi = static_cast<double>(n);
    f.y.fit(0, maxTotal > 0 ? maxTotal : 1, false, 0.12);
    f.yTicks = linearTicks(0, f.y.hi, 5);

    double lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
    bool has[3] = {false, false, false};
    for (const LineSeries& sr : overlays) {
        const int ax = sr.axis < 1 || sr.axis > 2 ? 1 : sr.axis;
        for (double v : sr.y) {
            if (!std::isfinite(v)) continue;
            if (!has[ax]) { lo[ax] = hi[ax] = v; has[ax] = true; }
            lo[ax] = std::min(lo[ax], v);
            hi[ax] = std::max(hi[ax], v);
        }
    }
    Axis* axes[3] = {&f.y, &f.y2, &f.y3};
    // Overlay axes track their own range instead of starting at zero: N50 and
    // unitig count move by a few percent across rounds, and a zero-based axis
    // would render both as flat lines.
    for (int i = 1; i < 3; ++i) {
        if (has[i]) axes[i]->fitRange(lo[i], hi[i], 0.22);
    }
    for (int i = 1; i < 3; ++i) {
        if (!has[i]) continue;
        std::vector<double>& dst = (i == 1) ? f.y2Ticks : f.y3Ticks;
        for (double t : f.yTicks) {
            const double g = frac(f.y.px(t) - f.y.p0, f.y.p1 - f.y.p0);
            dst.push_back(axes[i]->lo + g * (axes[i]->hi - axes[i]->lo));
        }
    }
    for (const LineSeries& sr : overlays) {
        if (sr.axis == 1) f.y2Color = sr.color;
        if (sr.axis == 2) f.y3Color = sr.color;
    }

    const double bw = (f.right() - f.left()) / static_cast<double>(n);
    const double base = f.y.px(0);
    for (size_t i = 0; i < n; ++i) {
        f.xTicks.push_back(static_cast<double>(i) + 0.5);
        f.xTickLabels.push_back(groupLabels[i]);
        double acc = 0;
        for (const StackSeries& sr : stack) {
            const double v = i < sr.v.size() ? sr.v[i] : 0.0;
            if (!(v > 0)) continue;
            const double y0 = f.y.px(acc);
            const double y1 = f.y.px(acc + v);
            const std::string tip = groupLabels[i] + " \xE2\x80\x94 " + sr.name + ": " +
                                    fmtInt(v) + " " + stackUnit;
            f.marks += svgRect(f.x.px(static_cast<double>(i)) + bw * 0.16, y1, bw * 0.68, y0 - y1,
                               "bar f" + std::to_string(sr.color), tip);
            acc += v;
        }
        if (totals[i] <= 0) {
            f.marks += svgRect(f.x.px(static_cast<double>(i)) + bw * 0.16, base - 1.5, bw * 0.68, 1.5,
                               "bar f6", groupLabels[i] + " \xE2\x80\x94 nothing removed");
        }
    }
    for (const StackSeries& sr : stack) f.legend += legendSwatch(sr.color, sr.name);

    for (const LineSeries& sr : overlays) {
        const int ax = sr.axis < 1 || sr.axis > 2 ? 1 : sr.axis;
        const Axis& ya = *axes[ax];
        std::vector<double> px, py;
        for (size_t i = 0; i < n && i < sr.y.size(); ++i) {
            if (!std::isfinite(sr.y[i])) continue;
            px.push_back(f.x.px(static_cast<double>(i) + 0.5));
            py.push_back(ya.px(sr.y[i]));
        }
        f.overlay += svgPolyline(px, py, "ln k" + std::to_string(sr.color));
        for (size_t i = 0; i < n && i < sr.y.size(); ++i) {
            if (!std::isfinite(sr.y[i])) continue;
            std::string c = "<circle" + at("cx", f.x.px(static_cast<double>(i) + 0.5)) +
                            at("cy", ya.px(sr.y[i])) + at("r", 3.0) +
                            at("class", "pt f" + std::to_string(sr.color));
            if (i < sr.tips.size() && !sr.tips[i].empty()) c += at("data-tip", sr.tips[i]);
            c += "/>";
            f.overlay += c;
        }
        f.legend += legendSwatch(sr.color, sr.name, true);
    }
    return f.render();
}

struct ScatterPoint {
    double x = 0, y = 0, r = 3;
    std::string tip;
};

std::string svgScatter(ChartOpts opts, const std::vector<ScatterPoint>& pts, bool logX, bool logY,
                       int color, const std::vector<double>& hGuides,
                       const std::vector<std::string>& hGuideLabels) {
    if (pts.empty()) return noDataChart(opts.title, "No contigs to plot.");
    double xlo = pts[0].x, xhi = pts[0].x, ylo = pts[0].y, yhi = pts[0].y;
    for (const ScatterPoint& p : pts) {
        xlo = std::min(xlo, p.x); xhi = std::max(xhi, p.x);
        ylo = std::min(ylo, p.y); yhi = std::max(yhi, p.y);
    }
    for (double g : hGuides) { ylo = std::min(ylo, g); yhi = std::max(yhi, g); }

    Frame f(std::move(opts));
    f.x.fit(logX ? xlo : std::min(0.0, xlo), xhi, logX, 0.06);
    f.y.fit(logY ? ylo : std::min(0.0, ylo), yhi, logY, 0.10);
    f.xTicks = f.x.log ? logTicks(f.x.lo, f.x.hi) : linearTicks(f.x.lo, f.x.hi, 6);
    f.yTicks = f.y.log ? logTicks(f.y.lo, f.y.hi) : linearTicks(f.y.lo, f.y.hi, 5);
    f.xGrid = true;

    for (size_t i = 0; i < hGuides.size(); ++i) {
        const double py = f.y.px(hGuides[i]);
        f.overlay += "<line class=\"vmark k6\"" + at("x1", f.left()) + at("y1", py) +
                     at("x2", f.right()) + at("y2", py) + "/>";
        if (i < hGuideLabels.size()) {
            f.overlay += "<text class=\"vlab t6\"" + at("x", f.right() - 4) + at("y", py - 4) +
                         " text-anchor=\"end\">" + hGuideLabels[i] + "</text>";
        }
    }
    for (const ScatterPoint& p : pts) {
        std::string c = "<circle" + at("cx", f.x.px(p.x)) + at("cy", f.y.px(p.y)) + at("r", p.r) +
                        at("class", "dot f" + std::to_string(color));
        if (!p.tip.empty()) c += at("data-tip", p.tip);
        c += "/>";
        f.marks += c;
    }
    return f.render();
}

// One horizontal stacked bar: the whole wall clock, split by stage.
struct TimeSeg {
    std::string name;
    double seconds = 0;
    int color = 1;
};

std::string svgTimeline(const std::vector<TimeSeg>& segs, double total) {
    double sum = 0;
    for (const TimeSeg& s : segs) sum += s.seconds;
    if (!(sum > 0)) {
        return noDataChart("Where the wall clock went",
                           "No per-stage timings were recorded for this run.");
    }
    ChartOpts o;
    o.title = "Where the wall clock went";
    o.subtitle = "Stacked to the full run time. Hover any block for its share.";
    o.xlab = "elapsed (seconds)";
    o.w = 820; o.h = 102;
    o.mL = 16; o.mR = 18; o.mT = 14; o.mB = 42;
    o.minWidth = 480;

    Frame f(o);
    f.x.lo = 0;
    f.x.hi = std::max(sum, total > 0 ? total : sum);
    f.xTicks = linearTicks(0, f.x.hi, 7);
    for (double t : f.xTicks) f.xTickLabels.push_back(fmtNum(t, t < 10 ? 1 : 0));

    const double barTop = 24, barH = 34;
    double acc = 0;
    for (const TimeSeg& s : segs) {
        if (!(s.seconds > 0)) continue;
        const double a = f.x.px(acc), b = f.x.px(acc + s.seconds);
        const std::string tip = s.name + " \xE2\x80\x94 " + fmtDuration(s.seconds) + " (" +
                                fmtPct(s.seconds, f.x.hi) + " of the run)";
        f.marks += svgRect(a, barTop, std::max(b - a, 0.7), barH, "bar f" + std::to_string(s.color), tip);
        acc += s.seconds;
    }
    f.marks += "<line class=\"axis\"" + at("x1", f.left()) + at("y1", barTop + barH) +
               at("x2", f.right()) + at("y2", barTop + barH) + "/>";
    for (double t : f.xTicks) {
        const double pxv = f.x.px(t);
        f.marks += "<line class=\"grid\"" + at("x1", pxv) + at("y1", barTop + barH) +
                   at("x2", pxv) + at("y2", barTop + barH + 4) + "/>";
    }
    for (size_t i = 0; i < f.xTicks.size(); ++i) {
        f.marks += "<text class=\"tick\"" + at("x", f.x.px(f.xTicks[i])) + at("y", barTop + barH + 16) +
                   " text-anchor=\"middle\">" + f.xTickLabels[i] + "</text>";
    }
    f.xTicks.clear();
    f.xTickLabels.clear();
    f.yTicks.clear();

    // The generic frame draws left/bottom axes that make no sense for a single
    // horizontal bar, so the bar draws its own and the frame's are suppressed.
    std::string svg = "<figure class=\"chart\"><figcaption class=\"chart-title\">" + o.title +
                      "</figcaption><p class=\"chart-sub\">" + o.subtitle + "</p>"
                      "<div class=\"chart-wrap\"><svg class=\"cv\" role=\"img\"" +
                      at("viewBox", "0 0 " + n2(o.w) + " " + n2(o.h)) +
                      at("style", "min-width:" + n2(o.minWidth) + "px") + "><title>" + o.title +
                      "</title>" + f.marks +
                      "<text class=\"axlab\"" + at("x", (f.left() + f.right()) / 2) +
                      at("y", o.h - 6) + " text-anchor=\"middle\">" + o.xlab + "</text></svg></div>";
    std::string legend;
    for (const TimeSeg& s : segs) {
        if (!(s.seconds > 0)) continue;
        legend += legendSwatch(s.color, htmlEscape(s.name) + " <b>" + fmtDuration(s.seconds) + "</b>");
    }
    svg += "<div class=\"legend\">" + legend + "</div></figure>";
    return svg;
}

// ------------------------------------------------------------ page fragments

std::string card(const std::string& label, const std::string& value, const std::string& sub,
                 const std::string& tip = std::string()) {
    std::string s = "<div class=\"card\"";
    if (!tip.empty()) s += at("data-tip", tip);
    s += "><div class=\"card-l\">" + label + "</div><div class=\"card-v\">" + value + "</div>";
    if (!sub.empty()) s += "<div class=\"card-s\">" + sub + "</div>";
    s += "</div>";
    return s;
}

std::string table(const std::vector<std::string>& headers,
                  const std::vector<std::vector<std::string>>& rows,
                  const std::string& cls = std::string()) {
    std::string s = "<div class=\"tablewrap\"><table";
    if (!cls.empty()) s += at("class", cls);
    s += "><thead><tr>";
    for (const std::string& h : headers) s += "<th>" + h + "</th>";
    s += "</tr></thead><tbody>";
    for (const std::vector<std::string>& r : rows) {
        s += "<tr>";
        for (const std::string& c : r) s += "<td>" + c + "</td>";
        s += "</tr>";
    }
    s += "</tbody></table></div>";
    return s;
}

std::string kv(const std::vector<std::pair<std::string, std::string>>& items) {
    std::string s = "<dl class=\"kv\">";
    for (const auto& p : items) s += "<dt>" + p.first + "</dt><dd>" + p.second + "</dd>";
    s += "</dl>";
    return s;
}

std::string notRun(const std::string& what, const std::string& why) {
    return "<div class=\"na\"><b>" + what + "</b> " + why + "</div>";
}

// ---------------------------------------------------------------- the verdict

struct Criterion {
    std::string name;
    std::string detail;
    double points = 0;
    double maxPoints = 0;
    std::string explain;
};

struct Verdict {
    std::string grade;
    std::string cls;
    double score = 0;
    std::vector<Criterion> criteria;
};

Verdict judge(const AssemblyReport& rep) {
    Verdict v;
    const double total = static_cast<double>(rep.totalLength);
    const double nContigs = static_cast<double>(rep.contigs.size());

    // Depth proxy: the coverage mode found at the first rung of the ladder is
    // the closest thing the run has to raw sequencing depth.
    double depth = rep.meanCoverage;
    std::string depthSource = "mean contig coverage";
    if (!rep.iterations.empty() && rep.iterations.front().peakCoverage > 0) {
        depth = rep.iterations.front().peakCoverage;
        depthSource = "k-mer coverage peak at k=" + std::to_string(rep.iterations.front().k);
    }

    {
        Criterion c;
        c.name = "Contiguity";
        c.maxPoints = 40;
        const double r = ratio(static_cast<double>(rep.n50), total);
        c.points = 40.0 * logSpan(r, 0.0005, 0.10);
        c.detail = "N50 " + fmtBp(static_cast<double>(rep.n50)) + " = " + fmtPct(static_cast<double>(rep.n50), total, 2) +
                   " of the assembly (L50 " + fmtInt(static_cast<double>(rep.l50)) + ")";
        c.explain = "Full marks when the N50 reaches 10% of the total length; none below 0.05%.";
        v.criteria.push_back(c);
    }
    {
        Criterion c;
        c.name = "Fragmentation";
        c.maxPoints = 20;
        c.points = 20.0 * (1.0 - logSpan(nContigs, 50, 2000));
        c.detail = fmtInt(nContigs) + " contigs, largest " + fmtBp(static_cast<double>(rep.largest));
        c.explain = "Full marks at 50 contigs or fewer; none at 2,000 or more.";
        v.criteria.push_back(c);
    }
    {
        Criterion c;
        c.name = "Coverage";
        c.maxPoints = 25;
        c.points = 25.0 * clamp01((depth - 5.0) / 25.0);
        c.detail = fmtNum(depth, 1) + "x (" + depthSource + ")";
        c.explain = "Full marks at 30x or above; none at 5x or below, where a de Bruijn graph "
                    "cannot separate errors from real sequence.";
        v.criteria.push_back(c);
    }
    {
        Criterion c;
        c.name = "Gap content";
        c.maxPoints = 15;
        const double g = ratio(static_cast<double>(rep.gapBases), total);
        c.points = 15.0 * clamp01(1.0 - g / 0.01);
        c.detail = fmtInt(static_cast<double>(rep.gapBases)) + " N bases (" +
                   fmtPct(static_cast<double>(rep.gapBases), total, 3) + " of the assembly)";
        c.explain = "Full marks with no scaffold gaps; none once 1% of the assembly is N.";
        v.criteria.push_back(c);
    }

    for (const Criterion& c : v.criteria) v.score += c.points;
    if (rep.contigs.empty()) v.score = 0;

    if (v.score >= 85)      { v.grade = "Excellent"; v.cls = "ok"; }
    else if (v.score >= 70) { v.grade = "Good";      v.cls = "ok"; }
    else if (v.score >= 50) { v.grade = "Fair";      v.cls = "warn"; }
    else                    { v.grade = "Poor";      v.cls = "bad"; }
    return v;
}

// ------------------------------------------------------------------ CSS / JS

const char* kDarkVars =
    "--bg:#0d1117;--panel:#151b23;--panel2:#1b222c;--fg:#e5e9f0;--muted:#93a0b4;"
    "--line:#2a323d;--grid:#232b36;--accent:#6ea8fe;--ok:#3fb950;--warn:#d6a531;--bad:#f0603a;"
    "--s1:#6ea8fe;--s2:#f0883e;--s3:#4bc99b;--s4:#c58af9;--s5:#ff7b9c;--s6:#8b98ac;--s7:#e3c34e;"
    "--shade:rgba(110,168,254,.13);--tipbg:#0b0f14;--tipfg:#e5e9f0;"
    "--shadow:0 1px 2px rgba(0,0,0,.5);";

std::string css() {
    std::string s;
    s += ":root{--bg:#f6f7f9;--panel:#ffffff;--panel2:#f0f2f6;--fg:#161a1f;--muted:#5b6572;"
         "--line:#dde1e8;--grid:#e8ebf0;--accent:#1f5fd0;--ok:#137a3d;--warn:#9a6b06;--bad:#bd2f22;"
         "--s1:#1f5fd0;--s2:#d9701e;--s3:#0f8f6b;--s4:#8b46c4;--s5:#c33a5e;--s6:#6c7887;--s7:#a8860a;"
         "--shade:rgba(31,95,208,.10);--tipbg:#161a1f;--tipfg:#f6f7f9;"
         "--shadow:0 1px 2px rgba(16,24,40,.07),0 1px 3px rgba(16,24,40,.05);}";
    s += "@media (prefers-color-scheme:dark){:root:not([data-theme=\"light\"]){";
    s += kDarkVars;
    s += "}}";
    s += ":root[data-theme=\"dark\"]{";
    s += kDarkVars;
    s += "}";

    s += "*{box-sizing:border-box}"
         "html{-webkit-text-size-adjust:100%}"
         "body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.55 -apple-system,"
         "BlinkMacSystemFont,'Segoe UI',Roboto,'Helvetica Neue',Arial,sans-serif;}"
         ".wrap{max-width:1180px;margin:0 auto;padding:0 18px 72px}"
         "h1{font-size:26px;margin:0 0 4px;letter-spacing:-.01em}"
         "h2{font-size:19px;margin:0 0 4px;letter-spacing:-.005em}"
         "h3{font-size:15px;margin:20px 0 6px}"
         "h4{font-size:13.5px;margin:14px 0 4px;color:var(--muted);text-transform:uppercase;"
         "letter-spacing:.06em}"
         "p{margin:0 0 10px}"
         "a{color:var(--accent)}"
         "code,pre{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:12.5px}"
         // File paths and command lines have no spaces to wrap at.
         "code{overflow-wrap:anywhere}"
         "ul{padding-left:20px;margin:0 0 10px}"
         "pre{background:var(--panel2);border:1px solid var(--line);border-radius:8px;padding:10px 12px;"
         "overflow-x:auto;white-space:pre-wrap;word-break:break-all;margin:0 0 10px}";

    s += ".masthead{background:var(--panel);border-bottom:1px solid var(--line);padding:18px 0 14px}"
         ".masthead .wrap{padding-bottom:0}"
         ".mhrow{display:flex;flex-wrap:wrap;gap:12px;align-items:flex-start;justify-content:space-between}"
         ".sub{color:var(--muted);font-size:13.5px;margin:0}"
         ".pill{display:inline-block;padding:1px 8px;border-radius:999px;background:var(--panel2);"
         "border:1px solid var(--line);font-size:12px;color:var(--muted);margin-right:6px}"
         "button{font:inherit;font-size:13px;color:var(--fg);background:var(--panel2);"
         "border:1px solid var(--line);border-radius:8px;padding:5px 11px;cursor:pointer}"
         "button:hover{border-color:var(--accent)}";

    s += "nav.toc{position:sticky;top:0;z-index:40;background:var(--panel);"
         "border-bottom:1px solid var(--line);overflow-x:auto;-webkit-overflow-scrolling:touch}"
         "nav.toc ol{list-style:none;display:flex;gap:2px;margin:0;padding:6px 18px;white-space:nowrap;"
         "max-width:1180px;margin:0 auto}"
         "nav.toc a{display:block;padding:4px 9px;border-radius:7px;font-size:12.5px;"
         "text-decoration:none;color:var(--muted)}"
         "nav.toc a:hover{background:var(--panel2);color:var(--fg)}"
         "nav.toc a.on{background:var(--panel2);color:var(--accent);font-weight:600}";

    s += "section{background:var(--panel);border:1px solid var(--line);border-radius:12px;"
         "padding:18px;margin:20px 0;box-shadow:var(--shadow);scroll-margin-top:56px}"
         ".lede{color:var(--muted);font-size:13.5px;max-width:82ch;margin:0 0 14px}"
         ".secnum{color:var(--muted);font-weight:600;font-size:12.5px;letter-spacing:.08em}";

    s += ".cards{display:grid;grid-template-columns:repeat(auto-fill,minmax(148px,1fr));gap:10px;margin:12px 0}"
         ".card{background:var(--panel2);border:1px solid var(--line);border-radius:10px;padding:10px 12px}"
         ".card-l{font-size:11.5px;text-transform:uppercase;letter-spacing:.05em;color:var(--muted)}"
         ".card-v{font-size:20px;font-weight:650;letter-spacing:-.01em;margin-top:2px;"
         "font-variant-numeric:tabular-nums;overflow-wrap:anywhere}"
         ".card-s{font-size:12px;color:var(--muted);margin-top:1px}";

    // Grid children default to min-width:auto, which lets a wide chart or table
    // push the 1fr track past the viewport. Zeroing it makes the child scroll
    // inside its own container instead of scrolling the page.
    s += ".grid2{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:16px}"
         ".grid2>*{min-width:0}"
         "@media (max-width:900px){.grid2{grid-template-columns:1fr}"
         ".cards{grid-template-columns:repeat(auto-fill,minmax(132px,1fr))}"
         "section{padding:14px}.wrap{padding:0 12px 56px}}";

    s += ".tablewrap{overflow-x:auto;-webkit-overflow-scrolling:touch;margin:8px 0 12px;"
         "border:1px solid var(--line);border-radius:10px}"
         "table{border-collapse:collapse;width:100%;font-size:13px}"
         "th,td{padding:6px 10px;text-align:right;white-space:nowrap;border-bottom:1px solid var(--line)}"
         "th:first-child,td:first-child{text-align:left}"
         "thead th{background:var(--panel2);color:var(--muted);font-weight:600;font-size:12px;"
         "position:sticky;top:0}"
         "tbody tr:last-child td{border-bottom:0}"
         "tbody tr:hover{background:var(--panel2)}"
         "td.num,th.num{font-variant-numeric:tabular-nums}"
         "table.sortable thead th{cursor:pointer;-webkit-user-select:none;user-select:none}"
         "table.sortable thead th:hover{color:var(--fg)}"
         "table.sortable thead th .ar{opacity:.45;font-size:10px}";

    s += "dl.kv{display:grid;grid-template-columns:minmax(0,auto) minmax(0,1fr);gap:2px 14px;"
         "margin:0 0 10px;font-size:13.5px}"
         "dl.kv dt{color:var(--muted);overflow-wrap:anywhere}"
         "dl.kv dd{margin:0;font-variant-numeric:tabular-nums;overflow-wrap:anywhere}";

    s += ".chart{margin:14px 0 6px;padding:0}"
         ".chart-title{font-size:13.5px;font-weight:650;margin:0 0 2px}"
         ".chart-sub{font-size:12.5px;color:var(--muted);margin:0 0 6px;max-width:84ch}"
         ".chart-wrap{overflow-x:auto;-webkit-overflow-scrolling:touch}"
         ".cv{width:100%;height:auto;display:block}"
         ".cv .grid{stroke:var(--grid);stroke-width:1}"
         ".cv .axis{stroke:var(--line);stroke-width:1}"
         ".cv .tick{fill:var(--muted);font-size:10.5px;font-variant-numeric:tabular-nums}"
         ".cv .axlab{fill:var(--muted);font-size:11.5px;font-weight:600}"
         ".cv .vlab{font-size:10.5px;font-weight:700}"
         ".cv .vmark{stroke-dasharray:4 3;stroke-width:1.4;fill:none}"
         ".cv .bar{shape-rendering:crispEdges}"
         ".cv .bar:hover{opacity:.72}"
         ".cv .hit{fill:transparent;stroke:none}"
         ".cv .hit:hover{fill:var(--shade)}"
         ".cv .shade{fill:var(--shade);stroke:none}"
         ".cv .ln{fill:none;stroke-width:2;stroke-linejoin:round;stroke-linecap:round}"
         ".cv .pt{stroke:var(--panel);stroke-width:1}"
         ".cv .dot{fill-opacity:.62;stroke:none}"
         ".cv .dot:hover{fill-opacity:1}"
         ".cv .area{fill-opacity:.34;stroke:none}"
         ".nodata{font-size:13px;color:var(--muted);background:var(--panel2);border:1px dashed var(--line);"
         "border-radius:10px;padding:14px}";

    for (int i = 1; i <= 7; ++i) {
        const std::string n = std::to_string(i);
        s += ".f" + n + "{fill:var(--s" + n + ")}";
        s += ".k" + n + "{stroke:var(--s" + n + ")}";
        s += ".a" + n + "{fill:var(--s" + n + ")}";
        s += ".t" + n + "{fill:var(--s" + n + ")}";
        s += ".bg" + n + "{background:var(--s" + n + ")}";
    }

    s += ".legend{display:flex;flex-wrap:wrap;gap:4px 14px;font-size:12px;color:var(--muted);"
         "margin:6px 0 2px}"
         ".lg{display:inline-flex;align-items:center;gap:5px}"
         ".sw{width:11px;height:11px;border-radius:3px;display:inline-block;flex:none}"
         ".sw.swline{height:3px;border-radius:2px;width:14px}";

    s += ".verdict{display:flex;flex-wrap:wrap;gap:16px;align-items:center;"
         "background:var(--panel2);border:1px solid var(--line);border-radius:12px;padding:14px 16px}"
         ".vgrade{font-size:28px;font-weight:750;letter-spacing:-.02em}"
         ".vgrade.ok{color:var(--ok)}.vgrade.warn{color:var(--warn)}.vgrade.bad{color:var(--bad)}"
         ".vscore{font-size:12.5px;color:var(--muted)}"
         ".vreasons{flex:1 1 340px;min-width:0}"
         ".bar-row{display:grid;grid-template-columns:112px 1fr 64px;gap:8px;align-items:center;"
         "font-size:12.5px;margin:3px 0}"
         ".bar-track{height:8px;border-radius:999px;background:var(--line);overflow:hidden}"
         ".bar-fill{height:100%;border-radius:999px}"
         ".bar-fill.ok{background:var(--ok)}.bar-fill.warn{background:var(--warn)}"
         ".bar-fill.bad{background:var(--bad)}"
         ".bar-val{text-align:right;color:var(--muted);font-variant-numeric:tabular-nums}"
         "@media (max-width:560px){.bar-row{grid-template-columns:92px 1fr 54px}}";

    s += ".na{background:var(--panel2);border:1px dashed var(--line);border-radius:10px;"
         "padding:12px 14px;font-size:13.5px;color:var(--muted);margin:8px 0}"
         ".na b{color:var(--fg)}"
         ".note{font-size:12.5px;color:var(--muted);max-width:84ch}"
         "details{background:var(--panel2);border:1px solid var(--line);border-radius:10px;"
         "padding:10px 14px;margin:10px 0}"
         "details>summary{cursor:pointer;font-weight:600;font-size:13.5px}"
         "details[open]>summary{margin-bottom:8px}"
         ".kpanel{border:1px solid var(--line);border-radius:12px;padding:14px;margin:14px 0;"
         "background:var(--panel)}"
         ".kpanel>h3{margin-top:0;display:flex;flex-wrap:wrap;gap:8px;align-items:baseline}"
         ".kname{font-size:17px;font-weight:700;font-variant-numeric:tabular-nums}"
         ".delta{font-size:12px;color:var(--muted)}"
         ".up{color:var(--ok)}.down{color:var(--bad)}"
         "footer{color:var(--muted);font-size:12.5px;text-align:center;padding:8px 0 0}";

    s += "#tip{position:fixed;z-index:99;display:none;pointer-events:none;background:var(--tipbg);"
         "color:var(--tipfg);border-radius:7px;padding:5px 9px;font-size:12px;max-width:320px;"
         "box-shadow:0 4px 16px rgba(0,0,0,.28);line-height:1.4}";
    return s;
}

const char* kScript =
    "(function(){"
    "var R=document.documentElement,K='tessera-theme';"
    "function ap(t){if(t==='light'||t==='dark'){R.setAttribute('data-theme',t);}else{R.removeAttribute('data-theme');}}"
    "var sv=null;try{sv=localStorage.getItem(K);}catch(e){}ap(sv);"
    "var b=document.getElementById('theme-btn');"
    "function lab(){var t=R.getAttribute('data-theme');return t==='dark'?'dark':(t==='light'?'light':'auto');}"
    "function upd(){if(b)b.textContent='Theme: '+lab();}upd();"
    "if(b)b.addEventListener('click',function(){var t=R.getAttribute('data-theme');"
    "var n=t==='dark'?'light':(t==='light'?'':'dark');ap(n);"
    "try{if(n){localStorage.setItem(K,n);}else{localStorage.removeItem(K);}}catch(e){}upd();});"
    "})();"

    "(function(){"
    "var tip=document.getElementById('tip');if(!tip)return;var vis=false;"
    "function place(e){tip.style.left='0px';tip.style.top='0px';"
    "var r=tip.getBoundingClientRect(),x=e.clientX+14,y=e.clientY+18;"
    "if(x+r.width>window.innerWidth-8)x=window.innerWidth-8-r.width;"
    "if(y+r.height>window.innerHeight-8)y=e.clientY-r.height-12;"
    "if(x<4)x=4;if(y<4)y=4;tip.style.left=x+'px';tip.style.top=y+'px';}"
    "document.addEventListener('mousemove',function(e){"
    "var el=e.target,t=null;"
    "while(el&&el.nodeType===1){if(el.getAttribute&&el.getAttribute('data-tip')){t=el.getAttribute('data-tip');break;}el=el.parentNode;}"
    "if(t){tip.textContent=t;if(!vis){tip.style.display='block';vis=true;}place(e);}"
    "else if(vis){tip.style.display='none';vis=false;}},{passive:true});"
    "document.addEventListener('mouseleave',function(){tip.style.display='none';vis=false;});"
    "})();"

    "(function(){"
    "var ts=document.querySelectorAll('table.sortable');"
    "Array.prototype.forEach.call(ts,function(tb){"
    "var hs=tb.querySelectorAll('thead th');"
    "Array.prototype.forEach.call(hs,function(th,ci){"
    "th.addEventListener('click',function(){"
    "var body=tb.tBodies[0];if(!body)return;"
    "var rows=Array.prototype.slice.call(body.rows);"
    "var dir=th.getAttribute('data-dir')==='asc'?-1:1;"
    "Array.prototype.forEach.call(hs,function(o){o.removeAttribute('data-dir');"
    "var m=o.querySelector('.ar');if(m)m.textContent='';});"
    "th.setAttribute('data-dir',dir===1?'asc':'desc');"
    "var mk=th.querySelector('.ar');if(mk)mk.textContent=dir===1?' \\u25B2':' \\u25BC';"
    "rows.sort(function(a,b){"
    "var ca=a.cells[ci],cb=b.cells[ci];if(!ca||!cb)return 0;"
    "var va=ca.getAttribute('data-v'),vb=cb.getAttribute('data-v');"
    "if(va!==null&&vb!==null){var na=parseFloat(va),nb=parseFloat(vb);"
    "if(na<nb)return -dir;if(na>nb)return dir;return 0;}"
    "return ca.textContent.localeCompare(cb.textContent)*dir;});"
    "var fr=document.createDocumentFragment();"
    "rows.forEach(function(r){fr.appendChild(r);});body.appendChild(fr);});});});"
    "})();"

    "(function(){"
    "if(!window.IntersectionObserver)return;"
    "var links={},secs=document.querySelectorAll('section[id]');"
    "Array.prototype.forEach.call(document.querySelectorAll('nav.toc a'),function(a){"
    "links[a.getAttribute('href').slice(1)]=a;});"
    "var io=new IntersectionObserver(function(es){"
    "es.forEach(function(en){var a=links[en.target.id];if(!a)return;"
    "if(en.isIntersecting){Array.prototype.forEach.call(document.querySelectorAll('nav.toc a.on'),"
    "function(o){o.className='';});a.className='on';}});},"
    "{rootMargin:'-56px 0px -70% 0px',threshold:0});"
    "Array.prototype.forEach.call(secs,function(s){io.observe(s);});"
    "})();";

// ------------------------------------------------------------------ sections

struct Section {
    const char* id;
    const char* label;
};

const Section kSections[] = {
    {"run",         "1 Run"},
    {"summary",     "2 Summary"},
    {"timeline",    "3 Timeline"},
    {"correction",  "4 Correction"},
    {"ladder",      "5 k ladder"},
    {"simplify",    "6 Simplification"},
    {"resolve",     "7 Repeat resolution"},
    {"polish",      "8 Polishing"},
    {"composition", "9 Composition"},
    {"outputs",     "10 Outputs"},
    {"config",      "11 Configuration"},
};

std::string sectionOpen(const char* id, const char* number, const std::string& title,
                        const std::string& lede) {
    return "<section id=\"" + std::string(id) + "\"><div class=\"secnum\">" + number +
           "</div><h2>" + title + "</h2><p class=\"lede\">" + lede + "</p>";
}

// ---- 1. header + verdict

std::string sectionRun(const AssemblyReport& rep, const Verdict& v) {
    std::string s = sectionOpen("run", "SECTION 1", "The run and its verdict",
        "tessera assembles short reads by building a de Bruijn graph over a ladder of k-mer "
        "sizes, cleaning that graph, walking through repeats with paired-end evidence and "
        "polishing the result against the reads. Everything below is what this particular run "
        "observed about itself, in the order it happened.");

    s += "<div class=\"verdict\">";
    s += "<div><div class=\"vgrade " + v.cls + "\">" + v.grade + "</div>"
         "<div class=\"vscore\">quality score " + fmtNum(v.score, 0) + " / 100</div></div>";
    s += "<div class=\"vreasons\">";
    // Bars keep the pipeline-facing order; the reasons below are re-sorted.
    for (const Criterion& c : v.criteria) {
        const double f = ratio(c.points, c.maxPoints);
        const char* cls = f >= 0.8 ? "ok" : (f >= 0.5 ? "warn" : "bad");
        s += "<div class=\"bar-row\"" +
             at("data-tip", htmlEscape(c.name + ": " + c.detail + ". " + c.explain)) + ">";
        s += "<div>" + htmlEscape(c.name) + "</div>";
        s += "<div class=\"bar-track\"><div class=\"bar-fill " + std::string(cls) + "\" style=\"width:" +
             n2(100.0 * f) + "%\"></div></div>";
        s += "<div class=\"bar-val\">" + fmtNum(c.points, 0) + "/" + fmtNum(c.maxPoints, 0) + "</div>";
        s += "</div>";
    }
    s += "</div></div>";

    std::vector<const Criterion*> order;
    for (const Criterion& c : v.criteria) order.push_back(&c);
    std::sort(order.begin(), order.end(), [](const Criterion* a, const Criterion* b) {
        return (a->maxPoints - a->points) > (b->maxPoints - b->points);
    });
    size_t lossy = 0;
    for (const Criterion* c : order) if (c->maxPoints - c->points > 0.5) ++lossy;
    s += "<h4>";
    s += lossy ? "What holds the score back" : "What the score rests on";
    s += "</h4><ul class=\"note\">";
    for (size_t i = 0; i < order.size() && i < 3; ++i) {
        const Criterion* c = order[i];
        const double lost = c->maxPoints - c->points;
        std::string verb;
        if (lost > c->maxPoints * 0.4)      verb = "is the main weakness";
        else if (lost > c->maxPoints * 0.1) verb = "costs a few points";
        else if (lost > 0.5)                verb = "is very slightly short of full marks";
        else                                verb = "is at full marks";
        s += "<li><b>" + htmlEscape(c->name) + "</b> " + verb + " &mdash; " +
             htmlEscape(c->detail) + ".</li>";
    }
    s += "</ul>";

    s += "<details><summary>How this verdict is computed</summary><p class=\"note\">"
         "The score is a reference-free heuristic worth 100 points, weighted toward the two "
         "things that make an assembly usable: how long the pieces are, and whether the data "
         "had the depth to justify them.</p>";
    std::vector<std::vector<std::string>> rows;
    for (const Criterion& c : v.criteria) {
        rows.push_back({htmlEscape(c.name), fmtNum(c.maxPoints, 0),
                        fmtNum(c.points, 1), htmlEscape(c.explain)});
    }
    s += table({"Criterion", "Weight", "Scored", "Rule"}, rows);
    s += "<p class=\"note\">Bands: 85+ excellent, 70+ good, 50+ fair, below 50 poor. "
         "No reference genome is involved, so this cannot detect a misassembly, contamination "
         "or a missing replicon &mdash; for those, read the coverage and GC scatters in "
         "section 9 and open the graph in Bandage.</p></details>";

    // ---- run facts
    s += "<div class=\"grid2\"><div>";
    s += "<h4>Execution</h4>";
    std::vector<std::pair<std::string, std::string>> exec;
    exec.push_back({"Tool", "tessera " + htmlEscape(rep.version)});
    exec.push_back({"Mode", htmlEscape(rep.mode)});
    exec.push_back({"Started", rep.startedAt.empty() ? std::string("not recorded")
                                                     : htmlEscape(rep.startedAt)});
    exec.push_back({"Wall clock", fmtDuration(rep.totalSeconds)});
    exec.push_back({"Peak memory", rep.peakMemoryBytes > 0 ? fmtBytes(rep.peakMemoryBytes)
                                                           : std::string("not available")});
    exec.push_back({"Threads", fmtInt(rep.threads)});
    s += kv(exec);
    s += "<h4>Command</h4>";
    if (rep.command.empty()) {
        s += "<p class=\"note\">The command line was not recorded in this run's report.</p>";
    } else {
        s += "<pre>" + htmlEscape(rep.command) + "</pre>";
    }
    s += "</div><div>";
    s += "<h4>Input</h4>";
    std::vector<std::pair<std::string, std::string>> in;
    in.push_back({"Reads", fmtInt(static_cast<double>(rep.reads))});
    in.push_back({"Bases", fmtInt(static_cast<double>(rep.inputBases)) + " (" +
                           fmtBp(static_cast<double>(rep.inputBases)) + ")"});
    in.push_back({"Max read length", fmtInt(rep.maxReadLength) + " bp"});
    in.push_back({"Mean read length",
                  rep.reads ? fmtNum(ratio(static_cast<double>(rep.inputBases),
                                           static_cast<double>(rep.reads)), 1) + " bp"
                            : std::string("n/a")});
    in.push_back({"Layout", rep.paired ? "paired-end" : "single-end / unpaired"});
    if (rep.totalLength > 0) {
        in.push_back({"Raw depth", fmtNum(ratio(static_cast<double>(rep.inputBases),
                                                static_cast<double>(rep.totalLength)), 1) +
                                   "x over the assembled length"});
    }
    s += kv(in);
    s += "<h4>Files</h4>";
    if (rep.inputFiles.empty()) {
        s += "<p class=\"note\">No input file names were recorded.</p>";
    } else {
        s += "<ul class=\"note\">";
        for (const std::string& fpath : rep.inputFiles) {
            s += "<li><code>" + htmlEscape(fpath) + "</code></li>";
        }
        s += "</ul>";
    }
    s += "</div></div></section>";
    return s;
}

// ---- 2. summary cards

std::string sectionSummary(const AssemblyReport& rep) {
    std::string s = sectionOpen("summary", "SECTION 2", "The assembly at a glance",
        "N50 is the length such that half the assembly sits in contigs at least that long; N90 "
        "does the same at 90%, and L50 is how many contigs it takes to reach half. Longer N50 "
        "with fewer contigs means a more contiguous assembly.");

    const double total = static_cast<double>(rep.totalLength);
    s += "<div class=\"cards\">";
    s += card("Contigs", fmtInt(static_cast<double>(rep.contigs.size())),
              "L50 " + fmtInt(static_cast<double>(rep.l50)),
              "Sequences written to contigs.fasta");
    s += card("Total length", fmtBp(total), fmtInt(total) + " bp",
              "Sum of every contig, gap bases included");
    s += card("Largest", fmtBp(static_cast<double>(rep.largest)),
              fmtPct(static_cast<double>(rep.largest), total) + " of total",
              "Longest single contig");
    s += card("N50", fmtBp(static_cast<double>(rep.n50)),
              fmtInt(static_cast<double>(rep.n50)) + " bp",
              "Half the assembly lies in contigs at least this long");
    s += card("N90", fmtBp(static_cast<double>(rep.n90)),
              fmtInt(static_cast<double>(rep.n90)) + " bp",
              "90% of the assembly lies in contigs at least this long");
    s += card("GC", fmtPctValue(rep.gcPercent, 2), "of called bases",
              "G+C share, gap bases excluded");
    s += card("Mean coverage", fmtNum(rep.meanCoverage, 1) + "x", "length-weighted",
              "Mean k-mer depth across contigs, weighted by contig length");
    s += card("Gap bases", fmtInt(static_cast<double>(rep.gapBases)),
              fmtPct(static_cast<double>(rep.gapBases), total, 3) + " of total",
              "N runs inserted by scaffolding");
    s += card("GFA segments", rep.gfaSegments ? fmtInt(static_cast<double>(rep.gfaSegments))
                                              : std::string("&mdash;"),
              rep.gfaSegments ? fmtInt(static_cast<double>(rep.gfaLinks)) + " links"
                              : std::string("graph not written"),
              "Nodes and edges in assembly_graph.gfa");
    s += card("Runtime", fmtDuration(rep.totalSeconds),
              fmtInt(rep.threads) + " threads", "Total wall clock for the run");
    s += "</div>";
    s += "</section>";
    return s;
}

// ---- 3. timeline

std::string sectionTimeline(const AssemblyReport& rep) {
    std::string s = sectionOpen("timeline", "SECTION 3", "Where the time went",
        "Each k on the ladder costs three things: counting its k-mers, building the compacted "
        "graph from the ones that survive the abundance cutoff, and simplifying that graph. "
        "The stages around the ladder &mdash; correction, repeat resolution, polishing &mdash; "
        "each make one more pass over the reads.");

    std::vector<TimeSeg> segs;
    double accounted = 0;
    if (rep.correctionRun && rep.correctionSeconds > 0) {
        segs.push_back({"error correction", rep.correctionSeconds, 7});
        accounted += rep.correctionSeconds;
    }
    for (const KIteration& it : rep.iterations) {
        const std::string kk = "k=" + std::to_string(it.k);
        segs.push_back({kk + " count", it.countSeconds, 1});
        segs.push_back({kk + " graph", it.graphSeconds, 2});
        segs.push_back({kk + " simplify", it.simplifySeconds, 3});
        accounted += it.countSeconds + it.graphSeconds + it.simplifySeconds;
    }
    if (rep.resolveRun && rep.resolveSeconds > 0) {
        segs.push_back({"repeat resolution", rep.resolveSeconds, 4});
        accounted += rep.resolveSeconds;
    }
    if (rep.polishRun && rep.polishSeconds > 0) {
        segs.push_back({"polishing", rep.polishSeconds, 5});
        accounted += rep.polishSeconds;
    }
    const double other = rep.totalSeconds - accounted;
    if (other > 0.001) segs.push_back({"reading input + writing output", other, 6});

    s += svgTimeline(segs, rep.totalSeconds);

    std::vector<std::vector<std::string>> rows;
    for (const TimeSeg& t : segs) {
        if (!(t.seconds > 0)) continue;
        rows.push_back({htmlEscape(t.name), fmtNum(t.seconds, 2) + " s",
                        fmtPct(t.seconds, rep.totalSeconds)});
    }
    rows.push_back({"<b>total</b>", "<b>" + fmtNum(rep.totalSeconds, 2) + " s</b>", "<b>100%</b>"});
    s += table({"Stage", "Seconds", "Share"}, rows);
    s += "</section>";
    return s;
}

// ---- 4. read error correction

std::string sectionCorrection(const AssemblyReport& rep) {
    std::string s = sectionOpen("correction", "SECTION 4", "Read error correction",
        "A sequencing error creates a short run of k-mers that appear nowhere else in the data. "
        "This stage counts k-mers at the smallest k on the ladder, treats everything above the "
        "abundance cutoff as trusted, then anchors each read on a stretch of trusted k-mers and "
        "extends outward: any base that breaks the run is replaced with the one that restores "
        "it. Doing this first shrinks the error cloud the later, larger-k counters have to hold, "
        "and rescues real sequence that would otherwise fall below the cutoff.");

    if (!rep.correctionRun) {
        s += notRun("Not run.",
                    "Read error correction was disabled for this run (<code>--no-correct</code>), "
                    "so the raw reads went straight into the k-mer ladder. Errors then have to be "
                    "removed by the abundance cutoff and by graph simplification instead.");
        s += "</section>";
        return s;
    }

    const CorrectionStats& c = rep.correction;
    const double examined = static_cast<double>(c.readsExamined);
    s += "<div class=\"cards\">";
    s += card("Anchor k", fmtInt(rep.correctionK), "smallest k on the ladder",
              "Correction uses the smallest k so trusted k-mers stay dense even in thin coverage");
    s += card("Reads examined", fmtInt(examined), "", "Every read the run loaded");
    s += card("Reads corrected", fmtInt(static_cast<double>(c.readsCorrected)),
              fmtPct(static_cast<double>(c.readsCorrected), examined) + " of reads",
              "Reads in which at least one base changed");
    s += card("Bases corrected", fmtInt(static_cast<double>(c.basesCorrected)),
              fmtPct(static_cast<double>(c.basesCorrected), static_cast<double>(rep.inputBases), 4) +
                  " of all bases",
              "Total substitutions applied");
    s += card("No trusted anchor", fmtInt(static_cast<double>(c.readsUncorrectable)),
              fmtPct(static_cast<double>(c.readsUncorrectable), examined, 3) + " of reads",
              "Reads with no trusted k-mer to start from; left untouched");
    const double perRead = ratio(static_cast<double>(c.basesCorrected),
                                 static_cast<double>(c.readsCorrected));
    s += card("Per corrected read", fmtNum(perRead, 2) + " bp", "mean edits",
              "Mean number of bases changed in the reads that were changed at all");
    s += "</div>";

    const double rate = ratio(static_cast<double>(c.basesCorrected),
                              static_cast<double>(rep.inputBases));
    s += "<p class=\"note\">Correction rate: <b>" + fmtNum(rate * 1e6, 1) +
         "</b> bases changed per million sequenced bases (" + fmtPct(rate, 1.0, 4) +
         "). A read with no trusted anchor is usually one that is short, low quality throughout, "
         "or from a region too thinly covered for its k-mers to pass the cutoff; those reads are "
         "left exactly as they were rather than guessed at.</p>";
    s += "</section>";
    return s;
}

// ---- 5. the multi-k ladder

std::string kmerHistogram(const KIteration& it) {
    // Index 0 is the count-zero bucket and is always empty; the chart starts at
    // abundance 1 where the error cloud lives.
    std::vector<double> vals;
    for (size_t i = 1; i < it.countHistogram.size(); ++i) {
        vals.push_back(static_cast<double>(it.countHistogram[i]));
    }
    ChartOpts o;
    o.title = "k-mer abundance spectrum at k = " + std::to_string(it.k);
    o.subtitle = "Distinct k-mers at each observed multiplicity, log scale. The left spike is the "
                 "error cloud &mdash; k-mers seen once or twice because a base was miscalled. The "
                 "broad hump is real genomic sequence, centred on the sequencing depth. Everything "
                 "left of the cutoff is discarded.";
    o.xlab = "k-mer abundance (times the k-mer was observed)";
    o.ylab = "distinct k-mers (log)";
    o.h = 320;
    o.mL = 62;

    std::vector<HistMarker> marks;
    if (it.cutoff > 0) {
        HistMarker m;
        m.x = static_cast<double>(it.cutoff);
        m.label = "cutoff " + std::to_string(it.cutoff);
        m.color = 5;
        m.labelDy = 12;
        m.tip = "Abundance cutoff " + std::to_string(it.cutoff) +
                ": k-mers seen fewer than this many times are dropped before the graph is built.";
        marks.push_back(m);
    }
    if (it.peakCoverage > 0) {
        HistMarker m;
        m.x = it.peakCoverage;
        m.label = "peak " + fmtNum(it.peakCoverage, 0) + "x";
        m.color = 3;
        m.labelDy = 26;
        m.tip = "Coverage peak at " + fmtNum(it.peakCoverage, 0) +
                "x: the most common multiplicity among genuine k-mers, i.e. the depth of "
                "single-copy sequence at this k.";
        marks.push_back(m);
    }
    const double shadeTo = it.cutoff > 0 ? static_cast<double>(it.cutoff) : 0.0;
    return svgHistogram(o, vals, 1.0, 1.0, true, 1, "x", "distinct k-mers", marks,
                        shadeTo > 0 ? 1.0 : 0.0, shadeTo,
                        "Discarded: k-mers below the abundance cutoff, almost all of them "
                        "sequencing errors.");
}

std::string sectionLadder(const AssemblyReport& rep) {
    std::string s = sectionOpen("ladder", "SECTION 5", "The multi-k ladder",
        "A de Bruijn graph breaks the reads into overlapping k-mers and joins them where they "
        "overlap by k&minus;1 bases. The choice of k is a trade-off with no single right answer: "
        "a small k keeps the graph connected through low-coverage regions but tangles every "
        "repeat shorter than k, while a large k walks straight through those repeats but falls "
        "apart wherever coverage dips or an error survives. tessera therefore runs the whole "
        "build at several k in increasing order. The contigs from each rung are fed into the next "
        "one as extra, heavily weighted evidence &mdash; the <i>carry-over</i> &mdash; which keeps "
        "sequence that only the smaller k could reach above the abundance cutoff at the larger k.");

    if (rep.iterations.empty()) {
        s += notRun("No iterations ran.",
                    "Every k on the ladder was skipped, which happens when the reads are shorter "
                    "than the smallest k requested.");
        s += "</section>";
        return s;
    }

    // Ladder-wide evolution.
    std::vector<double> xs;
    std::vector<std::string> xlabs;
    LineSeries n50s, units;
    n50s.name = "N50 after simplification (bp)"; n50s.color = 1; n50s.axis = 0;
    units.name = "unitigs after simplification"; units.color = 2; units.axis = 1;
    std::vector<std::vector<std::string>> ladderRows;
    for (size_t i = 0; i < rep.iterations.size(); ++i) {
        const KIteration& it = rep.iterations[i];
        xs.push_back(static_cast<double>(i));
        xlabs.push_back("k=" + std::to_string(it.k));
        n50s.y.push_back(static_cast<double>(it.n50Final));
        n50s.tips.push_back("k=" + std::to_string(it.k) + " \xE2\x80\x94 N50 " +
                            fmtInt(static_cast<double>(it.n50Final)) + " bp");
        units.y.push_back(static_cast<double>(it.unitigsFinal));
        units.tips.push_back("k=" + std::to_string(it.k) + " \xE2\x80\x94 " +
                             fmtInt(static_cast<double>(it.unitigsFinal)) + " unitigs");
        ladderRows.push_back({std::to_string(it.k),
                              fmtInt(it.cutoff),
                              fmtNum(it.peakCoverage, 0) + "x",
                              fmtInt(static_cast<double>(it.solidKmers)),
                              fmtInt(static_cast<double>(it.unitigsFinal)),
                              fmtInt(static_cast<double>(it.n50Final)),
                              fmtInt(static_cast<double>(it.lengthFinal)),
                              fmtInt(static_cast<double>(it.carryOverContigs))});
    }
    ChartOpts lo;
    lo.title = "How the graph changes as k grows";
    lo.subtitle = "N50 in bases on the left axis, unitig count on the right. A healthy ladder "
                  "climbs in N50 and falls in unitig count: the repeats that tangled the graph at "
                  "small k are being walked straight through at large k. Check the total length "
                  "column in the table below at the same time &mdash; if it collapses, the larger "
                  "k is losing sequence rather than resolving it.";
    lo.xlab = "k-mer size";
    lo.ylab = "N50 (bp)";
    lo.y2lab = "unitigs";
    lo.mR = 54;
    lo.h = 330;
    s += svgLineChart(lo, xs, {n50s, units}, xlabs);
    s += table({"k", "Cutoff", "Peak", "Solid k-mers", "Unitigs", "N50 (bp)", "Graph total (bp)",
                "Carry-over"}, ladderRows);

    // Per-k detail.
    for (const KIteration& it : rep.iterations) {
        s += "<div class=\"kpanel\"><h3><span class=\"kname\">k = " + std::to_string(it.k) +
             "</span>";
        s += "<span class=\"pill\">cutoff " + fmtInt(it.cutoff) + "</span>";
        s += "<span class=\"pill\">peak " + fmtNum(it.peakCoverage, 0) + "x</span>";
        if (it.carryOverContigs > 0) {
            s += "<span class=\"pill\">" + fmtInt(static_cast<double>(it.carryOverContigs)) +
                 " carry-over contigs</span>";
        } else {
            s += "<span class=\"pill\">first rung, no carry-over</span>";
        }
        s += "</h3>";

        s += kmerHistogram(it);

        s += "<div class=\"grid2\"><div><h4>k-mers</h4>";
        std::vector<std::pair<std::string, std::string>> a;
        a.push_back({"Instances counted", fmtInt(static_cast<double>(it.totalKmers))});
        a.push_back({"Distinct", fmtInt(static_cast<double>(it.distinctKmers))});
        a.push_back({"Solid (kept)", fmtInt(static_cast<double>(it.solidKmers)) + " &middot; " +
                                     fmtPct(static_cast<double>(it.solidKmers),
                                            static_cast<double>(it.distinctKmers))});
        a.push_back({"Dropped as error", fmtInt(static_cast<double>(it.distinctKmers) -
                                                static_cast<double>(it.solidKmers))});
        a.push_back({"Abundance cutoff", fmtInt(it.cutoff) + "x"});
        a.push_back({"Coverage peak", fmtNum(it.peakCoverage, 1) + "x"});
        a.push_back({"Median unitig coverage", fmtNum(it.medianCoverage, 1) + "x"});
        a.push_back({"Carry-over contigs", fmtInt(static_cast<double>(it.carryOverContigs))});
        s += kv(a);
        s += "</div><div><h4>Graph, before and after simplification</h4>";
        std::vector<std::vector<std::string>> rows;
        auto deltaCell = [](double before, double after, bool moreIsBetter) {
            const double d = after - before;
            const char* cls = (d > 0) == moreIsBetter ? "up" : "down";
            std::string txt = (d >= 0 ? "+" : "") + fmtInt(d);
            if (before > 0) txt += " (" + (d >= 0 ? std::string("+") : std::string()) +
                                   fmtNum(100.0 * d / before, 1) + "%)";
            if (d == 0) return std::string("<span class=\"delta\">no change</span>");
            return "<span class=\"delta " + std::string(cls) + "\">" + txt + "</span>";
        };
        rows.push_back({"Unitigs", fmtInt(static_cast<double>(it.unitigsBuilt)),
                        fmtInt(static_cast<double>(it.unitigsFinal)),
                        deltaCell(static_cast<double>(it.unitigsBuilt),
                                  static_cast<double>(it.unitigsFinal), false)});
        rows.push_back({"N50", fmtInt(static_cast<double>(it.n50Built)) + " bp",
                        fmtInt(static_cast<double>(it.n50Final)) + " bp",
                        deltaCell(static_cast<double>(it.n50Built),
                                  static_cast<double>(it.n50Final), true)});
        rows.push_back({"Total length", fmtInt(static_cast<double>(it.lengthBuilt)) + " bp",
                        fmtInt(static_cast<double>(it.lengthFinal)) + " bp",
                        deltaCell(static_cast<double>(it.lengthBuilt),
                                  static_cast<double>(it.lengthFinal), false)});
        s += table({"", "Built", "Simplified", "Change"}, rows);
        s += "<p class=\"note\">Time at this k: counting " + fmtNum(it.countSeconds, 2) +
             " s, graph " + fmtNum(it.graphSeconds, 2) + " s, simplification " +
             fmtNum(it.simplifySeconds, 2) + " s.</p>";
        s += "</div></div></div>";
    }

    s += "<p class=\"note\">Watch the abundance cutoff and coverage peak across the rungs. The "
         "peak falls as k rises &mdash; a read of length L contributes L&minus;k+1 k-mers, so "
         "each one is seen fewer times &mdash; and the cutoff follows it down. When the peak "
         "approaches the cutoff, the largest k is running out of signal and is unlikely to be "
         "adding contiguity.</p>";
    s += "</section>";
    return s;
}

// ---- 6. simplification

std::string sectionSimplify(const AssemblyReport& rep) {
    std::string s = sectionOpen("simplify", "SECTION 6", "Graph simplification",
        "The raw graph is full of structures that are artefacts rather than biology. Four "
        "operations run in rotation until nothing changes: <b>tip clipping</b> removes short "
        "dead-end branches whose coverage is weak next to the alternative at the branch point; "
        "<b>bubble popping</b> collapses two near-identical parallel paths, keeping the better "
        "supported one, but only when the losing side is clearly too shallow to be a real second "
        "repeat copy; <b>chimera removal</b> drops low-coverage unitigs that bridge two "
        "well-covered regions, the classic erroneous connection; and <b>isolated removal</b> "
        "deletes short, low-coverage fragments that connect to nothing. Each round is more "
        "aggressive than the last, so genuinely thin but real sequence gets a chance to be joined "
        "into something defensible before anything is willing to delete it.");

    bool any = false;
    for (const KIteration& it : rep.iterations) if (!it.rounds.empty()) any = true;
    if (!any) {
        s += notRun("No simplification rounds were recorded.",
                    "This happens when no k-mer iteration ran, or when the graph converged "
                    "before the first round completed.");
        s += "</section>";
        return s;
    }

    for (const KIteration& it : rep.iterations) {
        if (it.rounds.empty()) continue;
        std::vector<std::string> labels;
        StackSeries tips, bubbles, chimeras, isolated;
        tips.name = "tips"; tips.color = 1;
        bubbles.name = "bubbles"; bubbles.color = 2;
        chimeras.name = "chimeras"; chimeras.color = 5;
        isolated.name = "isolated"; isolated.color = 7;
        LineSeries n50Line, unitLine;
        n50Line.name = "N50 (bp)"; n50Line.color = 3; n50Line.axis = 1;
        unitLine.name = "unitigs"; unitLine.color = 6; unitLine.axis = 2;

        size_t totTips = 0, totBub = 0, totChi = 0, totIso = 0, totMerged = 0;
        std::vector<std::vector<std::string>> rows;
        for (const SimplifyRoundStats& r : it.rounds) {
            labels.push_back(std::to_string(r.round));
            tips.v.push_back(static_cast<double>(r.tipsRemoved));
            bubbles.v.push_back(static_cast<double>(r.bubblesPopped));
            chimeras.v.push_back(static_cast<double>(r.chimerasRemoved));
            isolated.v.push_back(static_cast<double>(r.isolatedRemoved));
            n50Line.y.push_back(static_cast<double>(r.n50));
            n50Line.tips.push_back("round " + std::to_string(r.round) + " \xE2\x80\x94 N50 " +
                                   fmtInt(static_cast<double>(r.n50)) + " bp");
            unitLine.y.push_back(static_cast<double>(r.unitigs));
            unitLine.tips.push_back("round " + std::to_string(r.round) + " \xE2\x80\x94 " +
                                    fmtInt(static_cast<double>(r.unitigs)) + " unitigs");
            totTips += r.tipsRemoved; totBub += r.bubblesPopped;
            totChi += r.chimerasRemoved; totIso += r.isolatedRemoved;
            totMerged += r.merged;
            rows.push_back({std::to_string(r.round),
                            fmtInt(static_cast<double>(r.tipsRemoved)),
                            fmtInt(static_cast<double>(r.bubblesPopped)),
                            fmtInt(static_cast<double>(r.chimerasRemoved)),
                            fmtInt(static_cast<double>(r.isolatedRemoved)),
                            fmtInt(static_cast<double>(r.merged)),
                            fmtInt(static_cast<double>(r.unitigs)),
                            fmtInt(static_cast<double>(r.n50)),
                            fmtInt(static_cast<double>(r.totalLength))});
        }
        rows.push_back({"<b>total</b>",
                        "<b>" + fmtInt(static_cast<double>(totTips)) + "</b>",
                        "<b>" + fmtInt(static_cast<double>(totBub)) + "</b>",
                        "<b>" + fmtInt(static_cast<double>(totChi)) + "</b>",
                        "<b>" + fmtInt(static_cast<double>(totIso)) + "</b>",
                        "<b>" + fmtInt(static_cast<double>(totMerged)) + "</b>", "", "", ""});

        s += "<div class=\"kpanel\"><h3><span class=\"kname\">k = " + std::to_string(it.k) +
             "</span><span class=\"pill\">" + std::to_string(it.rounds.size()) +
             " rounds</span></h3>";
        ChartOpts o;
        o.title = "Simplification rounds at k = " + std::to_string(it.k);
        o.subtitle = "Bars count graph elements removed in each round (left axis). The two lines "
                     "track what that bought: N50 in bases and the number of unitigs left, each "
                     "with its own scale in the matching colour on the right.";
        o.xlab = "simplification round";
        o.ylab = "elements removed";
        o.y2lab = "N50 (bp)";
        o.y3lab = "unitigs";
        o.mR = 92;
        o.h = 320;
        s += svgStackedBar(o, labels, {tips, bubbles, chimeras, isolated},
                           {n50Line, unitLine}, "removed");
        s += table({"Round", "Tips", "Bubbles", "Chimeras", "Isolated", "Merges", "Unitigs",
                    "N50 (bp)", "Total (bp)"}, rows);
        s += "</div>";
    }
    s += "<p class=\"note\">Merges are not deletions: after every removal the graph is "
         "re-compacted, and any chain of nodes that is now non-branching becomes a single longer "
         "unitig. That is where the N50 gain actually comes from.</p>";
    s += "</section>";
    return s;
}

// ---- 7. repeat resolution

std::string sectionResolve(const AssemblyReport& rep) {
    std::string s = sectionOpen("resolve", "SECTION 7", "Paired-end repeat resolution",
        "The compacted graph stops wherever a repeat creates a branch, even when the correct way "
        "through is obvious given the read pairs. This stage anchors reads onto unitigs using "
        "short probe k-mers, learns the fragment-length distribution from pairs that land on a "
        "single unitig, and then extends chains through the graph. Support is scored along whole "
        "candidate paths <i>through</i> the repeat rather than one branch at a time, and a pair "
        "only counts when the fragment length it would imply for that specific path is plausible "
        "&mdash; which is what separates a real continuation from a coincidental link. When two "
        "candidates are within the tie ratio of each other the extension is refused: forcing a "
        "near-tie is exactly how misassemblies get made.");

    if (!rep.resolveRun) {
        if (!rep.paired) {
            s += notRun("Not applicable.",
                        "The run had no paired reads, so there is no fragment-length information "
                        "to resolve repeats with. Contigs were taken straight from the simplified "
                        "graph, which means every repeat longer than the largest k remains a break "
                        "point.");
        } else {
            s += notRun("Not run.",
                        "Repeat resolution was disabled for this run "
                        "(<code>--no-resolve</code>). Contigs were taken straight from the "
                        "simplified graph, so each one ends at the first unresolved branch.");
        }
        s += "</section>";
        return s;
    }

    const ResolveStats& r = rep.resolve;
    const InsertModel& m = r.insert;

    s += "<div class=\"cards\">";
    s += card("Reads anchored", fmtInt(static_cast<double>(r.readsMapped)),
              fmtPct(static_cast<double>(r.readsMapped), static_cast<double>(rep.reads)) +
                  " of reads",
              "Reads placed at a unique position on a unitig");
    s += card("Linking pairs", fmtInt(static_cast<double>(r.pairsLinking)),
              "mates on different unitigs",
              "Pairs whose two reads anchored to different unitigs; these carry the joining "
              "information");
    s += card("Distinct links", fmtInt(static_cast<double>(r.distinctLinks)),
              "oriented unitig pairs", "Unique oriented unitig-to-unitig connections observed");
    s += card("Paths built", fmtInt(static_cast<double>(r.pathsBuilt)), "resolved walks",
              "Chains of unitigs emitted as single sequences");
    s += card("Unitigs joined", fmtInt(static_cast<double>(r.unitigsJoined)), "merged into paths",
              "Unitigs absorbed into a longer path");
    s += card("Scaffold joins", fmtInt(static_cast<double>(r.scaffoldJoins)),
              fmtInt(static_cast<double>(r.gapBases)) + " N bases",
              "Ends joined across a gap because paired support exists but no path through the "
              "graph does");
    s += "</div>";

    if (!m.usable) {
        s += notRun("The fragment-length model could not be estimated.",
                    "Too few pairs landed on a single unitig to fit a distribution, so branch "
                    "decisions fell back to raw link counts and graph topology only. That is "
                    "weaker evidence, and the resolution below should be read with that in mind.");
    } else {
        s += "<div class=\"grid2\"><div>";
        std::vector<double> vals;
        const size_t cap = std::min<size_t>(rep.insertHistogram.size(),
                                            static_cast<size_t>(m.maxPlausible > 0
                                                ? m.maxPlausible * 1.15 + 50 : 2000));
        for (size_t i = 0; i < cap; ++i) vals.push_back(static_cast<double>(rep.insertHistogram[i]));
        if (vals.empty()) {
            s += noDataChart("Fragment-length distribution",
                             "No fragment-length histogram was recorded.");
        } else {
            ChartOpts o = narrowChart();
            o.title = "Fragment-length distribution";
            o.subtitle = "Measured from pairs whose reads both land on the same unitig, where the "
                         "distance between them is known exactly. The shaded band is the window "
                         "of lengths the resolver will accept as plausible (mean &plusmn; 4 SD); "
                         "a pair implying a fragment outside it contributes nothing to a path's "
                         "score.";
            o.xlab = "fragment length (bp)";
            o.ylab = "read pairs";
            std::vector<HistMarker> marks;
            HistMarker mm;
            mm.x = m.mean;
            mm.label = "mean " + fmtNum(m.mean, 0) + " bp";
            mm.color = 5;
            mm.labelDy = 12;
            mm.tip = "Mean fragment length " + fmtNum(m.mean, 0) + " bp, SD " +
                     fmtNum(m.stddev, 0) + " bp, from " + fmtInt(static_cast<double>(m.observations)) +
                     " pairs.";
            marks.push_back(mm);
            s += svgHistogram(o, vals, 0.0, 1.0, false, 3, "bp", "pairs", marks,
                              static_cast<double>(m.minPlausible),
                              static_cast<double>(m.maxPlausible),
                              "Plausible window " + fmtInt(m.minPlausible) + "\xE2\x80\x93" +
                              fmtInt(m.maxPlausible) + " bp (mean \xC2\xB1 4 SD).");
        }
        s += "</div><div><h4>Fragment model</h4>";
        std::vector<std::pair<std::string, std::string>> im;
        im.push_back({"Mean", fmtNum(m.mean, 1) + " bp"});
        im.push_back({"Std deviation", fmtNum(m.stddev, 1) + " bp"});
        im.push_back({"Coefficient of variation", fmtPct(m.stddev, m.mean)});
        im.push_back({"Plausible window", fmtInt(m.minPlausible) + " &ndash; " +
                                          fmtInt(m.maxPlausible) + " bp"});
        im.push_back({"Observations", fmtInt(static_cast<double>(m.observations)) + " pairs"});
        s += kv(im);
        s += "<p class=\"note\">The window is the mean plus or minus four standard deviations. It "
             "is deliberately generous: its job is to reject links that would imply an absurd "
             "fragment, not to model the library precisely. A wide distribution (high coefficient "
             "of variation) makes the test less discriminating, so fewer repeats get resolved.</p>";
        s += "<h4>Yield</h4>";
        std::vector<std::pair<std::string, std::string>> yy;
        yy.push_back({"Anchor rate", fmtPct(static_cast<double>(r.readsMapped),
                                            static_cast<double>(rep.reads), 2)});
        yy.push_back({"Pairs per distinct link",
                      fmtNum(ratio(static_cast<double>(r.pairsLinking),
                                   static_cast<double>(r.distinctLinks)), 1)});
        yy.push_back({"Unitigs per path",
                      fmtNum(ratio(static_cast<double>(r.unitigsJoined),
                                   static_cast<double>(r.pathsBuilt)), 1)});
        yy.push_back({"Time", fmtDuration(rep.resolveSeconds)});
        s += kv(yy);
        s += "</div></div>";
    }

    if (r.scaffoldJoins > 0) {
        s += "<p class=\"note\"><b>Scaffolding.</b> " +
             fmtInt(static_cast<double>(r.scaffoldJoins)) +
             " joins were made across gaps totalling " + fmtInt(static_cast<double>(r.gapBases)) +
             " N bases. These are ends with strong paired support but no path through the graph "
             "connecting them &mdash; the sequence in between was never assembled, so a run of Ns "
             "sized from the fragment model stands in for it. The order and orientation are "
             "supported by the reads; the bases themselves are not.</p>";
    } else {
        s += "<p class=\"note\">No scaffold joins were made, so every contig is contiguous "
             "sequence with no N placeholders.</p>";
    }
    s += "</section>";
    return s;
}

// ---- 8. polishing

std::string sectionPolish(const AssemblyReport& rep) {
    std::string s = sectionOpen("polish", "SECTION 8", "Consensus polishing",
        "Resolution stitches unitigs together, and any base the graph got wrong survives into the "
        "contigs. Polishing re-anchors the reads onto the finished sequence and takes a "
        "per-position majority, changing a base only where the depth and the agreement among "
        "reads both clear a threshold. This is what drives the residual mismatch rate toward "
        "zero without inventing anything.");

    if (!rep.polishRun) {
        s += notRun("Not run.",
                    "Polishing was skipped for this run &mdash; either explicitly "
                    "(<code>--no-polish</code> or <code>--polish-passes 0</code>) or because the "
                    "chosen mode does not polish. The contigs are the graph consensus as it came "
                    "out of the previous stage.");
        s += "</section>";
        return s;
    }

    const PolishStats& p = rep.polish;
    const double total = static_cast<double>(rep.totalLength);
    s += "<div class=\"cards\">";
    s += card("Reads used", fmtInt(static_cast<double>(p.readsUsed)),
              fmtPct(static_cast<double>(p.readsUsed), static_cast<double>(rep.reads)) + " of reads",
              "Reads that anchored onto a contig and voted");
    s += card("Mean depth", fmtNum(p.meanDepth, 1) + "x", "read votes per base",
              "Mean number of reads covering a contig position");
    s += card("Positions covered", fmtInt(static_cast<double>(p.positionsCovered)),
              fmtPct(static_cast<double>(p.positionsCovered), total, 2) + " of the assembly",
              "Contig positions seen by at least one read");
    s += card("Low coverage", fmtInt(static_cast<double>(p.lowCoveragePositions)),
              fmtPct(static_cast<double>(p.lowCoveragePositions), total, 2) + " of the assembly",
              "Positions too thinly covered to be changed safely; left as they were");
    s += card("Bases changed", fmtInt(static_cast<double>(p.basesChanged)),
              fmtPct(static_cast<double>(p.basesChanged), total, 4) + " of the assembly",
              "Substitutions applied to the contigs");
    s += "</div>";

    if (p.basesChanged == 0) {
        s += "<p class=\"note\">No bases changed: the consensus coming out of the graph already "
             "agreed with the reads everywhere it was safe to check. That is the expected outcome "
             "on a clean isolate at good depth &mdash; it means the graph, not the polisher, did "
             "the work.</p>";
    } else {
        s += "<p class=\"note\">One change per " +
             fmtInt(ratio(total, static_cast<double>(p.basesChanged))) +
             " assembled bases. Each of these was a position where the graph consensus disagreed "
             "with a clear majority of the reads covering it.</p>";
    }
    const double uncovered = total - static_cast<double>(p.positionsCovered);
    if (uncovered > 0) {
        s += "<p class=\"note\">" + fmtInt(uncovered) + " positions (" +
             fmtPct(uncovered, total, 3) + ") were not covered by any anchored read and could not "
             "be checked at all. Scaffold gaps account for " +
             fmtInt(static_cast<double>(rep.gapBases)) + " of them.</p>";
    }
    s += "<p class=\"note\">Polishing took " + fmtDuration(rep.polishSeconds) + ".</p>";
    s += "</section>";
    return s;
}

// ---- 9. composition

std::string sectionComposition(const AssemblyReport& rep) {
    std::string s = sectionOpen("composition", "SECTION 9", "Assembly composition",
        "Four views of the finished contigs. Together they answer two questions that the summary "
        "numbers cannot: is the contiguity concentrated in a few long pieces or spread thin, and "
        "does everything in the assembly look like it came from the same genome.");

    if (rep.contigs.empty()) {
        s += notRun("No contigs were produced.",
                    "Nothing passed the minimum length filter, so there is nothing to plot.");
        s += "</section>";
        return s;
    }

    std::vector<size_t> lens;
    lens.reserve(rep.contigs.size());
    for (const ContigRecord& c : rep.contigs) lens.push_back(c.length);
    std::sort(lens.begin(), lens.end(), std::greater<size_t>());
    const double total = static_cast<double>(rep.totalLength);

    // --- length distribution, log-spaced bins
    {
        const double lmin = static_cast<double>(lens.back());
        const double lmax = static_cast<double>(lens.front());
        const double e0 = std::floor(std::log10(lmin > 0 ? lmin : 1) * 3.0) / 3.0;
        const double e1 = std::ceil(std::log10(lmax > 0 ? lmax : 1) * 3.0) / 3.0;
        std::vector<double> edges;
        for (double e = e0; e <= e1 + 1e-9 && edges.size() < 60; e += 1.0 / 3.0) {
            edges.push_back(std::pow(10.0, e));
        }
        if (edges.size() < 2) edges.push_back(edges.empty() ? 10.0 : edges[0] * 10);
        std::vector<double> counts(edges.size() - 1, 0.0);
        std::vector<double> bases(edges.size() - 1, 0.0);
        for (size_t L : lens) {
            for (size_t b = 0; b + 1 < edges.size(); ++b) {
                if (static_cast<double>(L) >= edges[b] &&
                    (static_cast<double>(L) < edges[b + 1] || b + 2 == edges.size())) {
                    counts[b] += 1;
                    bases[b] += static_cast<double>(L);
                    break;
                }
            }
        }
        std::vector<std::string> labels, tips;
        for (size_t b = 0; b + 1 < edges.size(); ++b) {
            labels.push_back(fmtTick(edges[b]));
            tips.push_back(fmtTick(edges[b]) + "\xE2\x80\x93" + fmtTick(edges[b + 1]) + " bp \xE2\x80\x94 " +
                           fmtInt(counts[b]) + " contigs, " + fmtBp(bases[b]) + " of sequence (" +
                           fmtPct(bases[b], total) + ")");
        }
        ChartOpts o = narrowChart();
        o.title = "Contig length distribution";
        o.subtitle = "Contigs per length bin, both axes logarithmic (bins are one third of a "
                     "decade wide). A tall left-hand column means many short fragments; check the "
                     "cumulative curve next to it to see whether they hold any real sequence.";
        o.xlab = "contig length (bp, log bins)";
        o.ylab = "contigs (log)";
        o.mB = 58;
        s += "<div class=\"grid2\"><div>";
        s += svgBarChart(o, labels, counts, true, 1, tips);
        s += "</div><div>";
    }

    // --- Nx cumulative curve
    {
        ChartOpts o = narrowChart();
        o.title = "Cumulative length curve (Nx)";
        o.subtitle = "Contigs laid out longest first: at x percent, the height is the length of "
                     "the contig that carries the assembly past that share of its total. N50 is "
                     "simply the height at 50%, N90 the height at 90% &mdash; a curve that stays "
                     "high for a long way is a contiguous assembly.";
        o.xlab = "percent of assembly (%)";
        o.ylab = "contig length (bp)";
        Frame f(o);
        f.x.lo = 0; f.x.hi = 100;
        f.y.fit(0, static_cast<double>(lens.front()), false, 0.08);
        f.xTicks = linearTicks(0, 100, 5);
        f.yTicks = linearTicks(0, f.y.hi, 5);
        f.xGrid = true;

        std::vector<double> px, py;
        double cum = 0;
        for (size_t i = 0; i < lens.size(); ++i) {
            const double x0 = 100.0 * ratio(cum, total);
            cum += static_cast<double>(lens[i]);
            const double x1 = 100.0 * ratio(cum, total);
            const double yv = f.y.px(static_cast<double>(lens[i]));
            px.push_back(f.x.px(x0)); py.push_back(yv);
            px.push_back(f.x.px(x1)); py.push_back(yv);
            if (lens.size() <= 400) {
                const std::string tip = "contig " + std::to_string(i + 1) + " \xE2\x80\x94 " +
                                        fmtInt(static_cast<double>(lens[i])) + " bp, spans " +
                                        fmtNum(x0, 1) + "\xE2\x80\x93" + fmtNum(x1, 1) +
                                        "% of the assembly";
                f.marks += svgRect(f.x.px(x0), f.top(), std::max(f.x.px(x1) - f.x.px(x0), 1.0),
                                   f.bottom() - f.top(), "hit", tip);
            }
        }
        f.marks += svgPolyline(px, py, "ln k1");

        auto guide = [&](double xPct, double yVal, const std::string& label, int color) {
            if (!(yVal > 0)) return;
            const double gx = f.x.px(xPct), gy = f.y.px(yVal);
            f.overlay += "<line class=\"vmark k" + std::to_string(color) + "\"" + at("x1", f.left()) +
                         at("y1", gy) + at("x2", gx) + at("y2", gy) + "/>";
            f.overlay += "<line class=\"vmark k" + std::to_string(color) + "\"" + at("x1", gx) +
                         at("y1", gy) + at("x2", gx) + at("y2", f.bottom()) + "/>";
            const bool flip = gx > f.left() + (f.right() - f.left()) * 0.68;
            f.overlay += "<text class=\"vlab t" + std::to_string(color) + "\"" +
                         at("x", gx + (flip ? -5 : 5)) + at("y", gy - 5) +
                         at("text-anchor", flip ? "end" : "start") + ">" + label + "</text>";
        };
        guide(50, static_cast<double>(rep.n50), "N50 " + fmtBp(static_cast<double>(rep.n50)), 5);
        guide(90, static_cast<double>(rep.n90), "N90 " + fmtBp(static_cast<double>(rep.n90)), 3);
        f.legend = legendSwatch(1, "contig length", true) + legendSwatch(5, "N50", true) +
                   legendSwatch(3, "N90", true);
        s += f.render();
        s += "</div></div>";
    }

    // --- scatters
    {
        const size_t cap = std::min<size_t>(rep.contigs.size(), 2000);
        std::vector<ScatterPoint> covLen, gcCov;
        double maxLen = 1;
        for (size_t i = 0; i < cap; ++i) maxLen = std::max(maxLen, static_cast<double>(rep.contigs[i].length));
        for (size_t i = 0; i < cap; ++i) {
            const ContigRecord& c = rep.contigs[i];
            const std::string name = "NODE_" + std::to_string(i + 1);
            const std::string tip = name + " \xE2\x80\x94 " + fmtInt(static_cast<double>(c.length)) +
                                    " bp, " + fmtNum(c.coverage, 1) + "x, GC " +
                                    fmtNum(c.gcPercent, 1) + "%" +
                                    (c.gapBases ? ", " + fmtInt(static_cast<double>(c.gapBases)) +
                                                  " N" : std::string());
            ScatterPoint p;
            p.x = static_cast<double>(c.length) > 0 ? static_cast<double>(c.length) : 1;
            p.y = c.coverage;
            p.r = 3.0;
            p.tip = tip;
            covLen.push_back(p);

            ScatterPoint q;
            q.x = c.gcPercent;
            q.y = c.coverage;
            // Area scales with contig length so the long contigs read as the ones
            // that actually define the assembly's centre of mass.
            q.r = 2.2 + 6.0 * std::sqrt(ratio(static_cast<double>(c.length), maxLen));
            q.tip = tip;
            gcCov.push_back(q);
        }

        s += "<div class=\"grid2\"><div>";
        {
            ChartOpts o = narrowChart();
            o.title = "Coverage against length";
            o.subtitle = "One point per contig. The long contigs should cluster on one horizontal "
                         "band &mdash; the chromosomal depth. A short contig sitting well above "
                         "that band is typically a plasmid or another multi-copy element; a "
                         "cluster well below it is usually contamination or a low-abundance "
                         "companion organism.";
            o.xlab = "contig length (bp, log)";
            o.ylab = "mean coverage (x)";
            s += svgScatter(o, covLen, true, false, 1, {rep.meanCoverage},
                            {"assembly mean " + fmtNum(rep.meanCoverage, 1) + "x"});
        }
        s += "</div><div>";
        {
            ChartOpts o = narrowChart();
            o.title = "GC against coverage";
            o.subtitle = "Same contigs, sized by length. One genome gives one tight cloud. A "
                         "second cloud at a different GC content is a strong sign of contamination "
                         "or of a plasmid with base composition unlike the chromosome.";
            o.xlab = "GC content (%)";
            o.ylab = "mean coverage (x)";
            s += svgScatter(o, gcCov, false, false, 4, {rep.meanCoverage},
                            {"assembly mean " + fmtNum(rep.meanCoverage, 1) + "x"});
        }
        s += "</div></div>";
        if (rep.contigs.size() > cap) {
            s += "<p class=\"note\">The scatters show the " + fmtInt(static_cast<double>(cap)) +
                 " longest contigs of " + fmtInt(static_cast<double>(rep.contigs.size())) +
                 " to keep the page small.</p>";
        }
    }

    // --- top 50 table
    {
        const size_t n = std::min<size_t>(rep.contigs.size(), 50);
        // Cells carry their own <td> here because each needs a numeric sort key
        // in data-v that the shared table() helper has no way to attach.
        auto cell = [](double v, const std::string& text) {
            return "<td class=\"num\" data-v=\"" + n2(v) + "\">" + text + "</td>";
        };
        std::string t = "<div class=\"tablewrap\"><table class=\"sortable\"><thead><tr>";
        const char* heads[] = {"Contig", "Length (bp)", "Coverage (x)", "GC (%)", "Gap bases",
                               "Cumulative (%)"};
        for (const char* h : heads) t += "<th>" + std::string(h) + "<span class=\"ar\"></span></th>";
        t += "</tr></thead><tbody>";
        double cum = 0;
        for (size_t i = 0; i < n; ++i) {
            const ContigRecord& c = rep.contigs[i];
            cum += static_cast<double>(c.length);
            t += "<tr><td data-v=\"" + std::to_string(i + 1) + "\">NODE_" + std::to_string(i + 1) +
                 "</td>";
            t += cell(static_cast<double>(c.length), fmtInt(static_cast<double>(c.length)));
            t += cell(c.coverage, fmtNum(c.coverage, 2));
            t += cell(c.gcPercent, fmtNum(c.gcPercent, 2));
            t += cell(static_cast<double>(c.gapBases), fmtInt(static_cast<double>(c.gapBases)));
            t += cell(100.0 * ratio(cum, total), fmtNum(100.0 * ratio(cum, total), 2));
            t += "</tr>";
        }
        t += "</tbody></table></div>";
        s += "<h4>Longest " + fmtInt(static_cast<double>(n)) + " contigs</h4>";
        s += "<p class=\"note\">Click any column header to sort. Names match the FASTA headers, "
             "which are <code>NODE_rank_length_L_cov_C</code>.</p>";
        s += t;
    }
    s += "</section>";
    return s;
}

// ---- 10. outputs

std::string sectionOutputs(const AssemblyReport& rep) {
    std::string s = sectionOpen("outputs", "SECTION 10", "How to use the outputs",
        "Four files land in the output directory. They describe the same assembly at different "
        "levels of commitment: the FASTA is the answer, the GFA is the evidence behind it.");

    s += "<h4>contigs.fasta</h4><p class=\"note\">The assembly, written longest first. Each header "
         "is <code>NODE_rank_length_L_cov_C</code>, where <code>L</code> is the sequence length and "
         "<code>C</code> the mean k-mer coverage &mdash; the same coverage plotted in section 9, so "
         "you can grep the headers to find candidate plasmids without leaving the shell. Runs of "
         "<code>N</code> are scaffold gaps: order and orientation are supported by the read pairs, "
         "the bases in between are not known.</p>";

    s += "<h4>assembly_graph.gfa</h4><p class=\"note\">The simplified unitig graph in GFA 1.0. "
         "<code>S</code> lines are unitigs (with <code>LN</code> length and <code>dp</code> depth "
         "tags), <code>L</code> lines are the (k&minus;1) overlaps between them, and <code>P</code> "
         "lines give each output contig as an ordered walk over those unitigs &mdash; so you can "
         "see exactly which graph nodes a contig traverses and where it passed through a repeat. "
         "Open it in Bandage to look at the topology: a contig that ends at a tangled knot of "
         "nodes is one whose repeat could not be resolved, and the <code>P</code> line tells you "
         "which nodes are involved.</p>";
    if (rep.gfaSegments == 0) {
        s += "<p class=\"note\">This run recorded no GFA segments, so the graph file was probably "
             "skipped with <code>--no-gfa</code>.</p>";
    } else {
        s += "<p class=\"note\">This run wrote " + fmtInt(static_cast<double>(rep.gfaSegments)) +
             " segments and " + fmtInt(static_cast<double>(rep.gfaLinks)) + " links for " +
             fmtInt(static_cast<double>(rep.contigs.size())) + " contigs. There are more segments "
             "than contigs because a contig is a walk over several unitigs, and because unitigs "
             "shorter than the minimum contig length stay in the graph without being emitted.</p>";
    }

    s += "<h4>report.json</h4><p class=\"note\">Every number on this page, machine-readable, with "
         "the same structure: run metadata, one object per k iteration including its full count "
         "histogram and per-round simplification stats, the resolution and polishing blocks, and a "
         "row per contig. Use it to compare runs or to plot the ladder yourself.</p>";

    s += "<h4>report.html</h4><p class=\"note\">This page. It is entirely self-contained &mdash; "
         "no scripts, styles, fonts or images are loaded from anywhere &mdash; so it can be "
         "archived or emailed and will still render years from now.</p>";

    s += "<p class=\"note\">A fifth file, <code>unitigs.fasta</code>, appears when the run was "
         "given <code>--unitigs</code>: the graph nodes before repeat resolution, useful when you "
         "want the raw evidence rather than the resolved walks.</p>";
    s += "</section>";
    return s;
}

// ---- 11. configuration

std::string sectionConfig(const AssemblyReport& rep) {
    std::string s = sectionOpen("config", "SECTION 11", "Effective configuration",
        "The settings this run actually used, as recorded in the report.");

    s += "<details open><summary>Settings and stages</summary>";
    std::vector<std::pair<std::string, std::string>> cfg;
    cfg.push_back({"Mode", htmlEscape(rep.mode)});
    cfg.push_back({"Threads", fmtInt(rep.threads)});
    std::string ladder;
    for (size_t i = 0; i < rep.iterations.size(); ++i) {
        if (i) ladder += ", ";
        ladder += std::to_string(rep.iterations[i].k);
    }
    cfg.push_back({"k ladder (as run)", ladder.empty() ? std::string("none") : ladder});
    cfg.push_back({"Read layout", rep.paired ? "paired-end" : "single-end"});
    cfg.push_back({"Error correction", rep.correctionRun
                       ? "on, at k=" + std::to_string(rep.correctionK)
                       : std::string("off")});
    cfg.push_back({"Repeat resolution", rep.resolveRun ? "on" : "off"});
    cfg.push_back({"Scaffolding", rep.resolveRun
                       ? (rep.resolve.scaffoldJoins > 0
                              ? "on, " + fmtInt(static_cast<double>(rep.resolve.scaffoldJoins)) +
                                    " joins made"
                              : std::string("on, no joins needed"))
                       : std::string("not applicable")});
    cfg.push_back({"Polishing", rep.polishRun ? "on" : "off"});
    cfg.push_back({"GFA output", rep.gfaSegments ? "written" : "not written"});
    s += kv(cfg);

    if (!rep.iterations.empty()) {
        std::vector<std::vector<std::string>> rows;
        for (const KIteration& it : rep.iterations) {
            rows.push_back({std::to_string(it.k), fmtInt(it.cutoff), fmtNum(it.peakCoverage, 1),
                            fmtInt(static_cast<double>(it.carryOverContigs)),
                            std::to_string(it.rounds.size())});
        }
        s += "<h4>Per-k thresholds chosen automatically</h4>";
        s += table({"k", "Abundance cutoff", "Coverage peak (x)", "Carry-over contigs",
                    "Simplify rounds"}, rows);
        s += "<p class=\"note\">The cutoff at each k is picked by finding the valley between the "
             "error spike and the coverage peak in that k's abundance spectrum, then capped at a "
             "quarter of the peak so a shallow run can never discard most of its real k-mers.</p>";
    }
    s += "<p class=\"note\">Tuning knobs that were left at their defaults are not recorded "
         "individually in the report; the mode above determines them, and any explicit override "
         "would appear on the command line in section 1.</p>";
    s += "</details></section>";
    return s;
}

// -------------------------------------------------------------------- assemble

std::string buildDocument(const AssemblyReport& rep) {
    const Verdict v = judge(rep);

    std::string s;
    s += "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">";
    s += "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">";
    s += "<title>tessera assembly report</title><style>";
    s += css();
    s += "</style></head><body><div id=\"tip\"></div>";

    s += "<header class=\"masthead\"><div class=\"wrap\"><div class=\"mhrow\"><div>";
    s += "<h1>tessera assembly report</h1>";
    s += "<p class=\"sub\"><span class=\"pill\">v" + htmlEscape(rep.version) + "</span>"
         "<span class=\"pill\">" + htmlEscape(rep.mode) + " mode</span>"
         "<span class=\"pill\">" + fmtInt(rep.threads) + " threads</span>"
         "<span class=\"pill\">" + fmtDuration(rep.totalSeconds) + "</span>";
    if (!rep.startedAt.empty()) s += "<span class=\"pill\">" + htmlEscape(rep.startedAt) + "</span>";
    s += "</p></div><div><button id=\"theme-btn\" type=\"button\">Theme</button></div>";
    s += "</div></div></header>";

    s += "<nav class=\"toc\" aria-label=\"Sections\"><ol>";
    for (const Section& sec : kSections) {
        s += "<li><a href=\"#" + std::string(sec.id) + "\">" + sec.label + "</a></li>";
    }
    s += "</ol></nav>";

    s += "<main class=\"wrap\">";
    s += sectionRun(rep, v);
    s += sectionSummary(rep);
    s += sectionTimeline(rep);
    s += sectionCorrection(rep);
    s += sectionLadder(rep);
    s += sectionSimplify(rep);
    s += sectionResolve(rep);
    s += sectionPolish(rep);
    s += sectionComposition(rep);
    s += sectionOutputs(rep);
    s += sectionConfig(rep);
    s += "<footer>Generated by tessera " + htmlEscape(rep.version) +
         " &middot; this page is self-contained and makes no network requests.</footer>";
    s += "</main><script>";
    s += kScript;
    s += "</script></body></html>\n";
    return s;
}

}  // namespace

bool writeHtmlReport(const std::string& path, const AssemblyReport& rep, std::string& error) {
    const std::string doc = buildDocument(rep);
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) { error = "cannot write " + path; return false; }
    const bool ok = std::fwrite(doc.data(), 1, doc.size(), f) == doc.size();
    std::fclose(f);
    if (!ok) { error = "write failed on " + path; return false; }
    return true;
}

}  // namespace ts
