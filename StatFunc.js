class StatFunc {
  static _lanczos_lgamma(x) {
    const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059,
               12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
    const g = 7;
    x -= 1;
    let sum = 0;
    for (let i = c.length - 1; i > 0; -- i) {
      sum += c[i] / (x + i);
    }
    sum += c[0];
    const base = x + g + 0.5;
    return (0.91893853320467274178 + Math.log(Math.abs(sum)) - base + Math.log(base) * (x + 0.5));
  }

  static _stirling_lgamma(x) {
    const c = [8.333333333333333333333333333333e-2, -2.77777777777777777777777777778e-3, 7.9365079365079365079365079365e-4,
               -5.9523809523809523809523809524e-4, 8.4175084175084175084175084175e-4, -1.91752691752691752691752691753e-3,
               6.41025641025641025641025641026e-3, -2.955065359477124183006535947712e-2, 1.7964437236883057316493849001589e-1,
               2.6440988521760198673491020321633e-1, 1.340286404416839199447895100069e+1];
    const eps = 1e-16;
    let sum = 0, x2 = x * x, d = x, i;
    for (i = 0; i < c.length; ++ i) {
      const v = c[i] / d;
      sum += v;
      if (Math.abs(v) < eps * sum) break;
      d *= x2;
    }
    return (0.91893853320467274178032973640562 + sum - x + Math.log(x) * (x - 0.5));
  } 

  static lgamma(x) {
    if (x < 0.5) {
      return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * x))) - StatFunc.lgamma(1 - x);
    } else if (x < 9) {
      return StatFunc._lanczos_lgamma(x);
    }
    return StatFunc._stirling_lgamma(x);
  }

  static gamma(x) {
    if (x < 0.5) return Math.PI / Math.sin(Math.PI * x) / StatFunc.gamma(1 - x);
    return Math.exp(StatFunc.lgamma(x));
  }

  static _beta_cdf_cf(x, a, b) {
    let ab = a + b, an = - ab * x / (a + 1), c = 1, d = 1 / (1 + an), f = d;
    const eps = 1e-16;
    let m;
    for (m = 1; m < 1000; ++ m) {
      const a2m = a + m + m;
      an = m * (b - m) * x / (a2m - 1) / a2m;
      c = 1 + an / c;
      d = 1 / (1 + an * d);
      f = c * d * f;
      an = - (a + m) * (ab + m) * x / a2m / (a2m + 1);
      c = 1 + an / c;
      d = 1 / (1 + an * d);
      const cd = c * d;
      f = cd * f;
      if (Math.abs(1 - cd) < eps) break;
    }
    return f;
  }

  static beta_cdf(x, a, b, right_tail = false) {
    let beta;
    if (x > (a + 1) / (a + b + 2)) {
      beta = StatFunc.beta_cdf(1 - x, b, a);
      right_tail = ! right_tail;
    } else {
      beta = Math.exp(StatFunc.lgamma(a + b) - StatFunc.lgamma(a) - StatFunc.lgamma(b) + a * Math.log(x) + b * Math.log(1 - x)) / a * StatFunc._beta_cdf_cf(x, a, b);
    }
    return (right_tail ? 1 - beta : beta);
  }

  static fdist_cdf(x, df1, df2, right_tail = false) {
    const df1x = df1 * x;
    return StatFunc.beta_cdf(df1x / (df1x + df2), df1 / 2, df2 / 2, right_tail);
  }

  static tdist_cdf(x, df, right_tail = false) {
    const p = 0.5 * StatFunc.beta_cdf(df / (x * x + df), df / 2, 0.5);
    return (right_tail ? p : 1 - p);
  }

  static normal_cdf(x, right_tail = false) {
    if (x < 0) {
      x = - x;
      right_tail = ! right_tail;
    }
    const x2 = x * x;
    const p = 0.39894228040143268 / (x + 2.92678600515804815) * (x2 + 8.42742300458043240 * x + 18.38871225773938487) /
              (x2 + 5.81582518933527391 * x + 8.97280659046817350) * (x2 + 7.30756258553673541 * x + 18.25323235347346525) /
              (x2 + 5.70347935898051437 * x + 10.27157061171363079) * (x2 + 5.66479518878470765 * x + 18.61193318971775795) /
              (x2 + 5.51862483025707963 * x + 12.72323261907760928) * (x2 + 4.91396098895240075 * x + 24.14804072812762821) /
              (x2 + 5.26184239579604207 * x + 16.88639562007936908) * (x2 + 3.83362947800146179 * x + 11.61511226260603247) /
              (x2 + 4.92081346632882033 * x + 24.12333774572479110) * Math.exp(- x2 / 2);
    return (right_tail ? p : 1 - p);
  }

  static normal_cdf_inv(p, right_tail = false) {
    const r1n = [-7.784894002430293e-3, -3.223964580411365e-1, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783];
    const r1d = [0, 7.784695709041462e-3, 3.224671290700398e-1, 2.445134137142996, 3.754408661907416, 1];
    const r2n = [-39.69683028665376, 220.9460984245205, -275.9285104469687, 138.3577518672690, -30.66479806614716, 2.506628277459239];
    const r2d = [-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155288572, 1];

    if (p > 0.97575) {
      p = 1 - p;
      right_tail = ! right_tail;
    }
    let x;
    if (p < 0.02425) {
      const z = Math.sqrt(-2 * Math.log(p));
      let a = r1n[0], b = r1d[0];
      for (let i = 1; i < r1n.length; ++ i) {
        a = r1n[i] + a * z;
        b = r1d[i] + b * z;
      }
      x = a / b;
    } else {
      const z = p - 0.5;
      const z2 = z * z;
      let a = r2n[0], b = r2d[0];
      for (let i = 1; i < r2n.length; ++ i) {
        a = r2n[i] + a * z2;
        b = r2d[i] + b * z2;
      }
      x = z * a / b;
    }
    return (right_tail ? - x : x);
  }

  static anova(data) {
    let SSB = 0, SSE = 0, k = data.length, n = 0, total_sum = 0, group_mean = [];
    if (k < 2) return undefined;
    for (const group of data) {
      n += group.length;
      const sum = group.reduce((acc, curr) => acc + curr, 0);
      total_sum += sum;
      const mean = sum / group.length;
      group_mean.push(mean);
      SSE = group.reduce((acc, curr) => (curr - mean) ** 2 + acc, SSE);
    }
    if (n <= k) return undefined;
    const total_mean = total_sum / n;
    for (let i = 0; i < k; ++ i) {
      SSB += (group_mean[i] - total_mean) ** 2 * data[i].length;
    }
    const df1 = k - 1;
    const df2 = n - k;
    const MSB = SSB / df1;
    const MSE = SSE / df2;
    const F = MSB / MSE;
    return [[SSB, df1, MSB, F, StatFunc.fdist_cdf(F, df1, df2, true)], [SSE, df2, MSE], [SSB + SSE, n - 1]];
  }
};

class StatGrid {
  constructor(table, cfg) {
    const data = cfg['data'];
    const vertical = ('vertical' in cfg ? cfg['vertical'] : true);
    const precision = ('precision' in cfg ? cfg['precision'] : undefined);
    let s = '';
    if ('caption' in cfg) {
      s += '<caption>' + cfg['caption'] + '</caption>';
    }
    const maxLength = data.reduce((acc, curr) => (curr.length > acc ? curr.length : acc), -Infinity);
    if ('colHeader' in cfg) {
      if ('rowHeader' in cfg && cfg['colHeader'].length <= (vertical ? data.length : maxLength)) s += "<th></th>";
      for (const hdr of cfg['colHeader']) {
        s += '<th>' + hdr + '</th>';
      }
    }
    if (vertical) {
      for (let j = 0; j < maxLength; ++ j) {
        s += '<tr>';
        if ('rowHeader' in cfg) {
          if (j < cfg['rowHeader'].length) {
            s += '<th style="text-align: left">' + cfg['rowHeader'][j] + '</th>';
          } else {
            s += '<th></th>';
          }
        }
        for (let i = 0; i < data.length; ++ i) {
          if (j >= data[i].length || isNaN(data[i][j])) {
            s += '<td></td>';
          } else {
            s += '<td>' + StatGrid.setPrecision(data[i][j], precision) + '</td>';
          }
        }
        s += '</tr>';
      }
    } else {
      for (let i = 0; i < data.length; ++ i) {
        s += '<tr>';
        if ('rowHeader' in cfg) {
          if (i < cfg['rowHeader'].length) {
            s += '<th style="text-align: left">' + cfg['rowHeader'][i] + '</th>';
          } else {
            s += '<th></th>';
          }
        }
        for (let j = 0; j < data[i].length; ++ j) {
          if (isNaN(data[i][j])) {
            s += '<td></td>';
          } else {
            s += '<td>' + StatGrid.setPrecision(data[i][j], precision) + '</td>';
          }
        }
        s += '</tr>';
      }
    }
    table.innerHTML = s;
  }

  static setPrecision(data, precision) {
    if (precision != undefined && typeof data == 'number' && ! Number.isInteger(data)) {
      return data.toPrecision(precision);
    }
    return data;
  }
}

class StatChart {
  constructor(ctx, cfg) {
    this.ctx = ctx;
    this.cfg = cfg;
    if (! ('options' in cfg)) {
      this.cfg['options'] = {};
    }
    this.width = ctx.canvas.width - 1;
    this.height = ctx.canvas.height - 1;
    ctx.canvas.style.width = ctx.canvas.width + 'px';
    ctx.canvas.style.height = ctx.canvas.height + 'px';
    // This changes the font size
    ctx.canvas.width *= window.devicePixelRatio;
    ctx.canvas.height *= window.devicePixelRatio;
    ctx.scale(window.devicePixelRatio, window.devicePixelRatio);
    ctx.font = "12px sans-serif";
    switch (this.cfg.type) {
      case 'line':
        this.lineplot();
        break;
      case 'boxplot':
        this.boxplot();
        break;
      case 'probability':
        this.probabilityplot();
        break;
    }
  }

  fillText(text, x, y) {
    this.ctx.fillText(text, Math.floor(x) + 0.5, Math.floor(y) + 0.5);
  }

  rect(x, y, width, height) {
    this.ctx.rect(Math.floor(x) + 0.5, Math.floor(y) + 0.5, Math.round(width), Math.round(height));
  }

  moveTo(x, y) {
    this.ctx.moveTo(Math.floor(x) + 0.5, Math.floor(y) + 0.5);
  }

  lineTo(x, y) {
    this.ctx.lineTo(Math.floor(x) + 0.5, Math.floor(y) + 0.5);
  }

  arc(x, y, r, startangle, endangle, anticlockwise = false) {
    this.ctx.arc(Math.floor(x) + 0.5, Math.floor(y) + 0.5, Math.round(r), startangle, endangle, anticlockwise);
  }

  static firstChoiceColors = [[0, 0, 255], [255, 0, 0], [0, 255, 0], [0, 255, 255], [255, 0, 255], [0, 128, 0],
                              [0, 0, 128], [255, 128, 0], [0, 128, 128], [0, 128, 255], [128, 0, 128], [128, 128, 0],
                              [128, 0, 255], [128, 0, 0], [190, 190, 0]];
  static getColor(index) {
    if (index < 15) return StatChart.firstChoiceColors[index];
    const h = 0.618033988749895 * index * 360 % 360;
    let v = 0.618033988749895 * index;
    v = 1 - 0.5 * (v - Math.trunc(v));
    return StatChart.hsv2rgb(h, 1, v);
  }

  static hsv2rgb(h, s, v) {
    const c = v * s;
    const x = c * (1 - Math.abs(h / 60 % 2 - 1));
    const m = v - c;
    let rgb;
    if (h < 60) {
      rgb = [c, x, 0];
    } else if (h < 120) {
      rgb = [x, c, 0];
    } else if (h < 180) {
      rgb = [0, c, x];
    } else if (h < 240) {
      rgb = [0, x, c];
    } else if (h < 300) {
      rgb = [x, 0, c];
    } else {
      rgb = [c, 0, x];
    }
    return [Math.round((rgb[0] + m) * 255), Math.round((rgb[1] + m) * 255), Math.round((rgb[2] + m) * 255)];
  }

  static getDataMinMax(datasets) {
    let min = Infinity, max = - Infinity; 
    for (const dataset of datasets) {
      if (dataset['data'].length > 0) {
        if (Array.isArray(dataset['data'][0])) {
          for (const data of dataset['data']) {
            [min, max] = data.reduce((acc, curr) => {if (curr < acc[0]) acc[0] = curr;
                                                     if (curr > acc[1]) acc[1] = curr; return acc;}, [min, max]);
          }
        } else {
          [min, max] = dataset['data'].reduce((acc, curr) => {if (curr < acc[0]) acc[0] = curr;
                                                              if (curr > acc[1]) acc[1] = curr; return acc;}, [min, max]);
        }
      }
    }
    return StatChart.adjustDataMinMax(min, max);
  }

  static adjustDataMinMax(min, max) {
    const adjustRatio = 0.05;
    const adjust = adjustRatio * (max - min);
    if (min >= 0) {
      max += adjust;
      min -= adjust;
      if (min < 0) {min = 0;}
    } else if (max <= 0) {
      max += adjust;
      min -= adjust;
      if (max > 0) {max = 0;}
    } else {
      max += adjust;
      min -= adjust;
    }
    return [min, max];
  }

  static getDataLabels(min, max) {
    const d = max - min;
    let powerOf10 = 1 - Math.floor(Math.log10(d));
    let n = Math.trunc(d * 10 ** powerOf10);
    const steps = [10, 20, 25, 5, 10];
    let i;
    if (n > 75) {
      i = 2;
    } else if (n > 60) {
      i = 1;
    } else if (n > 30) {
      i = 0;
    } else if (n > 15) {
      i = 3;
    } else {
      i = 2;
      ++ powerOf10;
    }
    let step = steps[i] / 10 ** powerOf10;
    let n1 = Math.ceil(min / step);
    let n2 = Math.floor(max / step);
    if (n2 - n1 > 5) {
      ++ i;
      if (i == 3) -- powerOf10;
      step = steps[i] * 0.1 ** powerOf10;
      n1 = Math.ceil(min / step);
      n2 = Math.floor(max / step);
    }
    let labels = new Array(n2 - n1 + 1);
    for (i = n1; i <= n2; ++ i) {
      labels[i - n1] = i * step;
    }
    return labels;
  }

  static shortFormat(data, precision = 5) {
    if (precision != undefined && typeof data == 'number' && ! Number.isInteger(data)) {
      let d = data.toPrecision(precision);
      let text = d.toString();
      if (/[.]/.test(text)) {
        if (/[eE]/.test(text)) {
        } else {
          text = text.replace(/[.]?0+$/, "");
        }
      }
      return text;
    }
    return data;
  }

  // Labels: {'texts': text labels, 'center': true/false, 'minmax': minimal and maximal values, 'values': value labels}
  drawLabelsAndGrids(xLabels, yLabels) {
    let [x1, y1, x2, y2] = [0, this.height, this.width, 0];

    const yMetrics = this.ctx.measureText("M");
    const textHeight = yMetrics.fontBoundingBoxAscent + yMetrics.fontBoundingBoxDescent;
    y1 -= 1.5 * textHeight;
    let labelList = ('texts' in yLabels ? yLabels['texts'] : yLabels['values']);
    let labelMaxWidth = labelList.reduce((acc, curr) => {const w = this.ctx.measureText(StatChart.shortFormat(curr)).width; return (w > acc ? w : acc);}, 0);
    x1 += labelMaxWidth;

    if ('title' in this.cfg['options'] && 'text' in this.cfg['options']['title']) {
      const metrics = this.ctx.measureText(this.cfg['options']['title']['text']);
      y2 += 5 + metrics.actualBoundingBoxAscent;
      this.ctx.beginPath()
      this.ctx.fillStyle = 'black';
      this.ctx.textAlign = 'center';
      this.ctx.textBaseline = 'alphabetic';
      this.fillText(this.cfg['options']['title']['text'], (x1 + x2) / 2, y2);
      y2 += 5 + metrics.actualBoundingBoxDescent;
    }

    let dx = 2;
    this.ctx.beginPath();
    this.ctx.fillStyle = 'black';
    this.ctx.strokeStyle = 'lightgray';
    this.ctx.textAlign = 'right';
    this.ctx.textBaseline = 'middle';
    for (let i = 0; i < labelList.length; ++ i) {
      const y = y1 + (yLabels['values'][i] - yLabels['minmax'][0]) / (yLabels['minmax'][1] - yLabels['minmax'][0]) * (y2 - y1);
      this.fillText(StatChart.shortFormat(labelList[i]), x1, y);
      this.moveTo(x1 + dx, y);
      this.lineTo(x2, y);
    }
    x1 += dx;
    labelList = ('texts' in xLabels ? xLabels['texts'] : xLabels['values']);
    this.ctx.textAlign = 'center';
    this.ctx.textBaseline = 'alphabetic';
    if ('values' in xLabels) {
      const y = y1 + 0.5 * textHeight + yMetrics.fontBoundingBoxAscent;
      for (let i = 0; i < labelList.length; ++ i) {
        const x = x1 + (xLabels['values'][i] - xLabels['minmax'][0]) / (xLabels['minmax'][1] - xLabels['minmax'][0]) * (x2 - x1);
        this.fillText(StatChart.shortFormat(labelList[i]), x, y);
        this.moveTo(x, y1);
        this.lineTo(x, y2);
      }
    } else {
      const center = ('center' in xLabels ? xLabels['center'] : true);
      labelMaxWidth = labelList.reduce((acc, curr) => {const w = this.ctx.measureText(StatChart.shortFormat(curr) + " ").width; return (w > acc ? w : acc);}, 0);
      const n = (center ? labelList.length : labelList.length - 1);
      const step = Math.ceil(labelMaxWidth * n / (x2 - x1));
      dx = (x2 - x1) / n;
      for (let i = 1; i < n; ++ i) {
        const x = x1 + i * dx;
        this.moveTo(x, y1);
        this.lineTo(x, y2);
      }
      const y = y1 + 0.5 * textHeight + yMetrics.fontBoundingBoxAscent;
      for (let i = 0, x = (center ? x1 + dx / 2 : x1); i < labelList.length; i += step) {
        this.fillText(StatChart.shortFormat(labelList[i]), x + i * dx, y);
        if (step > 1) {
          this.moveTo(x + i * dx, y1);
          this.lineTo(x + i * dx, y1 + 5);
        }
      }
    }
    this.ctx.stroke();
    return [x1, y1, x2, y2];
  }

  plotBox(data, x, dx, ymin, ymax, y1, y2) {
    if (data.length == 0) return;
    const sorted = data.toSorted((a, b) => a - b);
    let q1, q3, median;
    if (sorted.length < 3) {
      q1 = sorted[0];
      q3 = sorted[sorted.length - 1];
      median = (q1 + q3) / 2;
    } else {
      const n = sorted.length + 1;
      let v = n / 4;
      let k = Math.trunc(v);
      q1 = sorted[k - 1] + (v - k) * (sorted[k] - sorted[k - 1]);
      v = 3 * n / 4;
      k = Math.trunc(v);
      q3 = sorted[k - 1] + (v - k) * (sorted[k] - sorted[k - 1]);
      v = n / 2;
      k = Math.trunc(v);
      median = sorted[k - 1] + (v - k) * (sorted[k] - sorted[k - 1]);
    }
    const iqr15 = 1.5 * (q3 - q1);
    let minVal = q1 - iqr15, maxVal = q3 + iqr15;
    let count = 1;
    const ratio = (y2 - y1) / (ymax - ymin);
    for (let i = 0; i < sorted.length; ++ i) {
      if (sorted[i] > minVal) {
        minVal = sorted[i];
        break;
      } else if (sorted[i] == sorted[i + 1]) {
        ++ count;
      } else {
        count = 1;
      }
    }
    if (minVal > q1) minVal = q1;
    count = 1;
    for ( let i = sorted.length - 1; i >= 0; -- i) {
      if (sorted[i] < maxVal) {
        maxVal = sorted[i];
        break;
      } else if (sorted[i] == sorted[i - 1]) {
        ++ count;
      } else {
        const y = y1 + (sorted[i] - ymin) * ratio;
        let step = 6;
        if (step * (count - 1) > dx) step = dx / (count - 1);
        for (let j = 0, dj = (count - 1) / 2; j < count; ++ j) {
          const xc = x + (j - dj) * step;
          this.moveTo(xc + 2, y);
          this.arc(xc, y, 2, 0, 2 * Math.PI);
        }
        count = 1;
      }
    }
    if (maxVal < q3) maxVal = q3;
    let y = y1 + (median - ymin) * ratio;
    this.moveTo(x - dx / 2, y);
    this.lineTo(x + dx / 2, y);
    this.rect(x - dx / 2, y1 + (q3 - ymin) * ratio, dx, (q1 - q3) * ratio);
    y = y1 + (minVal - ymin) * ratio;
    this.moveTo(x, y1 + (q1 - ymin) * ratio);
    this.lineTo(x, y);
    this.moveTo(x - dx / 4, y);
    this.lineTo(x + dx / 4, y);
    y = y1 + (maxVal - ymin) * ratio;
    this.moveTo(x, y1 + (q3 - ymin) * ratio);
    this.lineTo(x, y);
    this.moveTo(x - dx / 4, y);
    this.lineTo(x + dx / 4, y);
  }

  boxplot() {
    const [ymin, ymax] = StatChart.getDataMinMax(this.cfg['data']['datasets']);
    if (ymin == Infinity) return;  // No point to draw
    const ylabels = StatChart.getDataLabels(ymin, ymax);
    const [x1, y1, x2, y2] = this.drawLabelsAndGrids({'texts': this.cfg['data']['labels']}, {'minmax': [ymin, ymax], 'values': ylabels});

    const dx = (x2 - x1) / this.cfg['data']['labels'].length;
    for (const dataset of this.cfg['data']['datasets']) {
      for (let i = 0; i < dataset['data'].length; ++ i) {
        this.ctx.beginPath();
        const [r, g, b] = StatChart.getColor(i);
        this.ctx.strokeStyle = `rgb(${r}, ${g}, ${b})`;
        this.plotBox(dataset['data'][i], x1 + (i + 0.5) * dx, 0.8 * dx, ymin, ymax, y1, y2);
        this.ctx.stroke();
      }
    }

    this.ctx.beginPath();
    this.rect(x1, y2, x2 - x1, y1 - y2);
    this.ctx.strokeStyle = 'black';
    this.ctx.stroke();
  }

  plotProb(data, xmin, xmax, ymin, ymax, x1, x2, y1, y2) {
    if (data.length == 0) return;
    const sorted = data.toSorted((a, b) => a - b);
    const xRatio = (x2 - x1) / (xmax - xmin), yRatio = (y2 - y1) / (ymax - ymin);
    const n1 = data.length + 1;
    this.moveTo(x1 + (sorted[0] - xmin) * xRatio, y1 + (StatFunc.normal_cdf_inv(1 / n1) - ymin) * yRatio);
    for (let i = 1; i < data.length; ++ i) {
      this.lineTo(x1 + (sorted[i] - xmin) * xRatio, y1 + (StatFunc.normal_cdf_inv((i + 1) / n1) - ymin) * yRatio);
    }
  }

  probabilityplot() {
    const [xmin, xmax] = StatChart.getDataMinMax(this.cfg['data']['datasets']);
    if (xmin == Infinity) return;
    const xlabels = StatChart.getDataLabels(xmin, xmax);
    let ymin = Infinity;
    for (const dataset of this.cfg['data']['datasets']) {
      ymin = dataset['data'].reduce((acc, curr) => {const v = 100 / (curr.length + 1); return (v < acc ? v : acc)}, ymin);
    }
    const yInv = StatFunc.normal_cdf_inv(ymin / 100);
    ymin = StatFunc.normal_cdf(StatChart.adjustDataMinMax(yInv, -yInv)[0]) * 100;
    let ylabels = [0, 25, 50, 75, 100];
    if (ymin <= 1) {
      ylabels[0] = 1;
      ylabels[4] = 99;
    } else if (ymin <= 2) {
      ylabels[0] = 2;
      ylabels[4] = 98;
    } else if (ymin <= 5) {
      ylabels[0] = 5;
      ylabels[4] = 95;
    } else if (ymin <=10) {
      ylabels[0] = 10;
      ylabels[4] = 90;
    } else {
      ylabels[0] = 25;
      ylabels[4] = 75;
    }
    ymin = StatFunc.normal_cdf_inv(ymin / 100);
    const ymax = - ymin;
    const [x1, y1, x2, y2] = this.drawLabelsAndGrids({'minmax': [xmin, xmax], 'values': xlabels}, {'texts': ylabels, 'minmax': [ymin, ymax], 'values': ylabels.map((item) => StatFunc.normal_cdf_inv(item / 100))});

    for (const dataset of this.cfg['data']['datasets']) {
      for (let i = 0; i < dataset['data'].length; ++ i) {
        this.ctx.beginPath();
        const [r, g, b] = StatChart.getColor(i);
        this.ctx.strokeStyle = `rgb(${r}, ${g}, ${b})`;
        this.plotProb(dataset['data'][i], xmin, xmax, ymin, ymax, x1, x2, y1, y2);
        this.ctx.stroke();
      }
    }

    this.ctx.beginPath();
    this.rect(x1, y2, x2 - x1, y1 - y2);
    this.ctx.strokeStyle = 'black';
    this.ctx.stroke();
  }

  plotLine(data, ymin, ymax, x1, x2, y1, y2) {
    if (data.length == 0) return;
    const dx = (x2 - x1) / (data.length - 1), yRatio = (y2 - y1) / (ymax - ymin);
    this.moveTo(x1, y1 + (data[0] - ymin) * yRatio);
    for (let i = 1; i < data.length; ++ i) {
      this.lineTo(x1 + i * dx, y1 + (data[i] - ymin) * yRatio);
    }
  }

  lineplot() {
    const [ymin, ymax] = StatChart.getDataMinMax(this.cfg['data']['datasets']);
    if (ymin == Infinity) return;  // No point to draw
    const ylabels = StatChart.getDataLabels(ymin, ymax);
    const [x1, y1, x2, y2] = this.drawLabelsAndGrids({'texts': this.cfg['data']['labels'], 'center': false}, {'minmax': [ymin, ymax], 'values': ylabels});

    for (const dataset of this.cfg['data']['datasets']) {
      for (let i = 0; i < dataset['data'].length; ++ i) {
        this.ctx.beginPath();
        const [r, g, b] = StatChart.getColor(i);
        this.ctx.strokeStyle = `rgb(${r}, ${g}, ${b})`;
        this.plotLine(dataset['data'][i], ymin, ymax, x1, x2, y1, y2);
        this.ctx.stroke();
      }
    }

    this.ctx.beginPath();
    this.rect(x1, y2, x2 - x1, y1 - y2);
    this.ctx.strokeStyle = 'black';
    this.ctx.stroke();
  }
};
