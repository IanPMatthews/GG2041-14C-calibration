let intcalCurve = [];
let calibrationChart = null;
let agedepthChart = null;
let lastCalibration = null;
let lastPDF = null;
let lastAgeDepthModel = [];

// Load IntCal20 JSON
fetch("simplified_intcal20.json")
  .then(r => r.json())
  .then(data => {
    intcalCurve = data.map(row => ({
      calendar_age: row.calBP,
      c14_age: row.C14Age,
      c14_sig: row.C14Sigma
    }));
    document.getElementById("curve-status").textContent =
      `Loaded IntCal20: ${intcalCurve.length} points`;
  })
  .catch(err => {
    document.getElementById("curve-status").textContent =
      `Error loading curve: ${err}`;
  });

// Linear interpolation helper
function linearInterp(x, x1, y1, x2, y2) {
  if (x2 === x1) return y1;
  return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}

// Interpolate IntCal in calendar-age space
function interpolateIntCalByCalendarAge(t) {
  if (t <= intcalCurve[0].calendar_age) return intcalCurve[0];
  if (t >= intcalCurve[intcalCurve.length - 1].calendar_age)
    return intcalCurve[intcalCurve.length - 1];

  for (let i = 0; i < intcalCurve.length - 1; i++) {
    const p1 = intcalCurve[i];
    const p2 = intcalCurve[i + 1];
    if (t >= p1.calendar_age && t <= p2.calendar_age) {
      return {
        c14_age: linearInterp(t, p1.calendar_age, p1.c14_age,
                                 p2.calendar_age, p2.c14_age),
        c14_sig: linearInterp(t, p1.calendar_age, p1.c14_sig,
                                 p2.calendar_age, p2.c14_sig)
      };
    }
  }
}

// Real calibrated PDF
function computeCalibratedPDF(rcAge, rcSigma) {
  const tMin = intcalCurve[0].calendar_age;
  const tMax = intcalCurve[intcalCurve.length - 1].calendar_age;

  const xs = [];
  const ys = [];

  for (let t = tMin; t <= tMax; t += 1) {
    const { c14_age, c14_sig } = interpolateIntCalByCalendarAge(t);
    const diff = rcAge - c14_age;
    const varTotal = rcSigma**2 + c14_sig**2;
    const pdf = Math.exp(-0.5 * diff*diff / varTotal) /
                Math.sqrt(2 * Math.PI * varTotal);
    xs.push(t);
    ys.push(pdf);
  }

  const sum = ys.reduce((a, b) => a + b, 0);
  return { xs, ys: ys.map(v => v / sum) };
}

// Posterior mean + sigma
function calibrateRadiocarbon(rcAge, rcSigma) {
  const pdf = computeCalibratedPDF(rcAge, rcSigma);
  const xs = pdf.xs;
  const ys = pdf.ys;

  const mean = xs.reduce((acc, x, i) => acc + x * ys[i], 0);
  const variance = xs.reduce((acc, x, i) => {
    const d = x - mean;
    return acc + d*d * ys[i];
  }, 0);

  lastPDF = pdf;
  return { calendarAge: mean, calendarSigmaApprox: Math.sqrt(variance) };
}
// Highest Posterior Density interval
function computeHPD(xs, ys, prob) {
  const idx = xs.map((_, i) => i).sort((a, b) => ys[b] - ys[a]);
  const included = new Array(xs.length).fill(false);
  let acc = 0;

  for (const i of idx) {
    acc += ys[i];
    included[i] = true;
    if (acc >= prob) break;
  }

  let low = null, high = null;
  for (let i = 0; i < xs.length; i++) {
    if (included[i]) {
      if (low === null) low = xs[i];
      high = xs[i];
    }
  }
  return { low, high };
}

// Summary statistics (mean, median, mode, HPDs)
function computeSummary(xs, ys) {
  const mean = xs.reduce((acc, x, i) => acc + x * ys[i], 0);

  // Median
  let cum = 0;
  let median = xs[xs.length - 1];
  for (let i = 0; i < xs.length; i++) {
    cum += ys[i];
    if (cum >= 0.5) {
      median = xs[i];
      break;
    }
  }

  // Mode
  let mode = xs[0];
  let maxy = ys[0];
  for (let i = 1; i < xs.length; i++) {
    if (ys[i] > maxy) {
      maxy = ys[i];
      mode = xs[i];
    }
  }

  return {
    mean,
    median,
    mode,
    hpd68: computeHPD(xs, ys, 0.682),
    hpd95: computeHPD(xs, ys, 0.954)
  };
}

// Render summary table
function renderSummary(summary) {
  const div = document.getElementById("summary-output");
  div.innerHTML = `
    <table>
      <tr><th>Statistic</th><th>Value (cal yr BP)</th></tr>
      <tr><td>Mean</td><td>${summary.mean.toFixed(1)}</td></tr>
      <tr><td>Median</td><td>${summary.median.toFixed(1)}</td></tr>
      <tr><td>Mode</td><td>${summary.mode.toFixed(1)}</td></tr>
      <tr><td>68.2% HPD</td><td>${summary.hpd68.low.toFixed(1)} – ${summary.hpd68.high.toFixed(1)}</td></tr>
      <tr><td>95.4% HPD</td><td>${summary.hpd95.low.toFixed(1)} – ${summary.hpd95.high.toFixed(1)}</td></tr>
    </table>
  `;
}
// plot cal pdf
function plotCalibrationPDF(xs, ys) {
  const ctx = document.getElementById("calibration-chart").getContext("2d");

  if (calibrationChart) calibrationChart.destroy();

  const summary = computeSummary(xs, ys);

// Density-based trimming
const maxPDF = Math.max(...ys);
let trimmed = xs.map((x, i) => ({ x, y: ys[i] }))
  .filter(p => p.y > maxPDF * 0.01);

// Extract trimmed x-range
let tx = trimmed.map(p => p.x);
let ty = trimmed.map(p => p.y);

// Add ±200 yr padding
const pad = 200;
const paddedLow = tx[0] - pad;
const paddedHigh = tx[tx.length - 1] + pad;

// Re-trim xs/ys to padded range
trimmed = xs.map((x, i) => ({ x, y: ys[i] }))
  .filter(p => p.x >= paddedLow && p.x <= paddedHigh);

tx = trimmed.map(p => p.x);
ty = trimmed.map(p => p.y);


  // HPD shading band
  const hpdBand = tx.map((x, i) => {
    return (x >= summary.hpd95.low && x <= summary.hpd95.high)
      ? ty[i]
      : 0;
  });

  calibrationChart = new Chart(ctx, {
    type: "line",
    data: {
      labels: tx,
      datasets: [
        {
          label: "95% HPD",
          data: hpdBand,
          borderWidth: 0,
          pointRadius: 0,
          backgroundColor: "rgba(0, 120, 212, 0.15)",
          fill: true
        },
        {
          label: "Calibrated PDF",
          data: ty,
          borderColor: "#005a9e",
          backgroundColor: "rgba(0, 90, 158, 0.10)",
          pointRadius: 0,
          borderWidth: 2
        }
      ]
    },
    options: {
      responsive: true,
      scales: {
        x: {
          title: { display: true, text: "Calendar age (yr BP)" },
          ticks: { color: "#333" }
        },
        y: {
          title: { display: true, text: "Probability density" },
          ticks: { color: "#333" }
        }
      },
      plugins: {
        legend: {
          labels: { color: "#333" }
        }
      }
    }
  });
}



// Parse tie points
function parseTiePoints(text) {
  return text.split("\n")
    .map(l => l.trim())
    .filter(l => l.length)
    .map(l => {
      const [d, a] = l.split(",").map(Number);
      return { depth: d, age: a };
    })
    .sort((a, b) => a.depth - b.depth);
}

// Build piecewise linear age–depth model
function buildAgeDepthModel(points, start, end, step) {
  const out = [];
  for (let d = start; d <= end; d += step) {
    let seg = null;
    for (let i = 0; i < points.length - 1; i++) {
      const p1 = points[i], p2 = points[i + 1];
      if (d >= p1.depth && d <= p2.depth) {
        seg = { p1, p2 };
        break;
      }
    }
    if (!seg) {
      out.push({ depth: d, age: null });
    } else {
      out.push({
        depth: d,
        age: linearInterp(d, seg.p1.depth, seg.p1.age, seg.p2.depth, seg.p2.age)
      });
    }
  }
  return out;
}

// Plot age–depth model
function plotAgeDepth(model) {
  const ctx = document.getElementById("agedepth-chart").getContext("2d");

  const depths = model.filter(r => r.age !== null).map(r => r.depth);
  const ages = model.filter(r => r.age !== null).map(r => r.age);

  if (agedepthChart) agedepthChart.destroy();

  agedepthChart = new Chart(ctx, {
    type: "line",
    data: {
      labels: depths,
      datasets: [{
        label: "Age–Depth",
        data: ages,
        borderColor: "#9bd4ff",
        backgroundColor: "rgba(155, 212, 255, 0.2)",
        pointRadius: 0,
        borderWidth: 2
      }]
    },
    options: {
      responsive: true,
      scales: {
        x: {
          title: { display: true, text: "Depth (cm)" },
          ticks: { color: "#ccc" }
        },
        y: {
          title: { display: true, text: "Age (yr BP)" },
          ticks: { color: "#ccc" }
        }
      },
      plugins: {
        legend: {
          labels: { color: "#ddd" }
        }
      }
    }
  });
}

// CSV download helper
function downloadCSV(filename, header, rows) {
  const csv = [header.join(",")]
    .concat(rows.map(r => r.map(v => String(v)).join(",")))
    .join("\n");

  const blob = new Blob([csv], { type: "text/csv" });
  const url = URL.createObjectURL(blob);

  const a = document.createElement("a");
  a.href = url;
  a.download = filename;

  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);

  URL.revokeObjectURL(url);
}

// UI wiring
document.addEventListener("DOMContentLoaded", () => {

  // Calibration button
  document.getElementById("calibrate-btn").onclick = () => {
    const rcAge = Number(document.getElementById("rc-age").value);
    const rcSigma = Number(document.getElementById("rc-sigma").value);
    const out = document.getElementById("calibration-output");

    if (!Number.isFinite(rcAge) || !Number.isFinite(rcSigma)) {
      out.textContent = "Enter numeric radiocarbon age and sigma.";
      return;
    }

    const pdf = computeCalibratedPDF(rcAge, rcSigma);
    const res = calibrateRadiocarbon(rcAge, rcSigma);

    lastCalibration = { rcAge, rcSigma, ...res };
    lastPDF = pdf;

    out.textContent =
      `Input: ${rcAge} ± ${rcSigma} 14C yr BP\n` +
      `Posterior mean: ${res.calendarAge.toFixed(1)} cal yr BP\n` +
      `Posterior σ: ${res.calendarSigmaApprox.toFixed(1)} yr`;

    const summary = computeSummary(pdf.xs, pdf.ys);
    renderSummary(summary);
    plotCalibrationPDF(pdf.xs, pdf.ys);
  };

  // Calibration CSV download
  document.getElementById("download-calibration-csv").onclick = () => {
    if (!lastCalibration) return;

    const rows = [[
      lastCalibration.rcAge,
      lastCalibration.rcSigma,
      lastCalibration.calendarAge,
      lastCalibration.calendarSigmaApprox
    ]];

    downloadCSV(
      "calibration_result.csv",
      ["rc_age_bp", "rc_sigma", "posterior_mean_cal_bp", "posterior_sigma"],
      rows
    );
  };

  // Age–depth model button
  document.getElementById("run-model-btn").onclick = () => {
    const points = parseTiePoints(document.getElementById("tiepoints").value);
    const start = Number(document.getElementById("model-start").value);
    const end = Number(document.getElementById("model-end").value);
    const step = Number(document.getElementById("model-step").value);

    const out = document.getElementById("agedepth-output");

    try {
      const model = buildAgeDepthModel(points, start, end, step);
      lastAgeDepthModel = model;

      out.textContent = "depth_cm, age_yr_BP\n" +
        model.map(r => `${r.depth}, ${r.age ?? "NA"}`).join("\n");

      plotAgeDepth(model);

    } catch (err) {
      out.textContent = `Error: ${err.message}`;
    }
  };

  // Age–depth CSV download
  document.getElementById("download-agedepth-csv").onclick = () => {
    if (!lastAgeDepthModel.length) return;

    const rows = lastAgeDepthModel.map(r => [r.depth, r.age ?? "NA"]);

    downloadCSV(
      "age_depth_model.csv",
      ["depth_cm", "age_yr_BP"],
      rows
    );
  };

});
