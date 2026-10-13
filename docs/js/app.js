/**
 * QGAT — OBLO Interactive Studio & Platform Logic
 * Gene Annotation, QTLdb Annotation, and Combined Analysis
 */

(function () {
  'use strict';

  // ---- Datasets & Presets ----
  const ANIMAL_PRESETS = {
    cattle: {
      name: "Cattle (Bos taurus)",
      assembly: "ARS-UCD1.2",
      range: "14\t1750000\t1850000\n6\t37900000\t38200000\n2\t6400000\t6500000\n20\t31800000\t32200000\n26\t21100000\t21200000",
      pos: "14\t1801116\n6\t38024921\n2\t6445100\n20\t31920050\n26\t21150000"
    },
    goat: {
      name: "Goat (Capra hircus)",
      assembly: "CHIR_1.0",
      range: "6\t86000000\t87500000\n3\t15000000\t16500000\n11\t42000000\t43500000",
      pos: "6\t87145000\n3\t15800000\n11\t42900000"
    },
    sheep: {
      name: "Sheep (Ovis aries)",
      assembly: "OAR3.1",
      range: "X\t48000000\t49500000\n6\t35000000\t36500000\n2\t118000000\t120000000",
      pos: "X\t48750000\n6\t35800000\n2\t119000000"
    },
    pig: {
      name: "Pig (Sus scrofa)",
      assembly: "Sscrofa11.1",
      range: "2\t1400000\t1600000\n6\t49000000\t51000000\n4\t62000000\t64000000",
      pos: "2\t1500000\n6\t50200000\n4\t63100000"
    },
    chicken: {
      name: "Chicken (Gallus gallus)",
      assembly: "GRCg4",
      range: "4\t72000000\t74000000\n1\t54000000\t56000000\n3\t31000000\t33000000",
      pos: "4\t73000000\n1\t55100000\n3\t32000000"
    },
    horse: {
      name: "Horse (Equus caballus)",
      assembly: "EquCab2",
      range: "18\t66000000\t67000000\n3\t104000000\t106000000\n9\t45000000\t47000000",
      pos: "18\t66500000\n3\t105000000\n9\t46200000"
    }
  };

  const LOCUS_PRESETS = {
    bovine_fat: {
      mode: 'gene',
      species: 'cattle',
      coords: "14\t1790000\t1820000\n26\t21140000\t21170000\n6\t37900000\t38100000"
    },
    goat_fleece: {
      mode: 'qtl',
      species: 'goat',
      coords: "6\t86000000\t87500000\n3\t15000000\t16500000"
    },
    sheep_ovulation: {
      mode: 'qtl',
      species: 'sheep',
      coords: "X\t48000000\t49500000\n6\t35000000\t36500000"
    },
    porcine_muscle: {
      mode: 'combined',
      species: 'pig',
      coords: "2\t1400000\t1600000\n6\t49000000\t51000000"
    }
  };

  // State
  const state = {
    currentMode: 'gene',
    selectedAnimal: 'cattle',
    geneFormat: 'range',
    qtlFormat: 'range',
    combinedFormat: 'range',
    activeResultsTab: 'genes',
    hasResults: false,
    resultsGenes: [],
    resultsQTLs: []
  };

  // Initialize
  function init() {
    setupNavbar();
    setupTypeSelector();
    setupFormatToggles();
    setupWindowChips();
    setupAnimalGrid();
    setupUploadZones();
    setupPresetButtons();
    setupExecutionButtons();
    setupResultsTabs();
    setupSearchFilters();
    setupDownloadButtons();
    setupCliTabs();
    setupModal();
    setupExampleDatasetLoader();

    // Ensure results start empty per user instruction
    setResultsState(false);
  }

  // Navbar behavior
  function setupNavbar() {
    const nav = document.getElementById('mainNav');
    window.addEventListener('scroll', () => {
      if (window.scrollY > 40) {
        nav.classList.add('scrolled');
      } else {
        nav.classList.remove('scrolled');
      }
    });

    const mobileBtn = document.getElementById('mobileMenuBtn');
    const navLinks = document.getElementById('navLinks');
    if (mobileBtn && navLinks) {
      mobileBtn.addEventListener('click', () => {
        navLinks.classList.toggle('open');
      });
    }
  }

  // 1. Type Selector (Gene / QTL / Combined)
  function setupTypeSelector() {
    const cards = document.querySelectorAll('.type-card');
    cards.forEach(card => {
      card.addEventListener('click', () => {
        cards.forEach(c => c.classList.remove('selected'));
        card.classList.add('selected');

        const mode = card.dataset.type;
        state.currentMode = mode;

        document.querySelectorAll('.input-panel').forEach(p => p.classList.remove('active'));

        if (mode === 'gene') {
          document.getElementById('panelGene').classList.add('active');
        } else if (mode === 'qtl') {
          document.getElementById('panelQTL').classList.add('active');
        } else if (mode === 'combined') {
          document.getElementById('panelCombined').classList.add('active');
        }
      });
    });
  }

  // 2. Format Toggles (Range vs Position-Only)
  function setupFormatToggles() {
    function bindToggle(toggleId, windowGroupId, textareaId, labelId, stateKey) {
      const toggle = document.getElementById(toggleId);
      const windowGroup = document.getElementById(windowGroupId);
      const textarea = document.getElementById(textareaId);
      const label = document.getElementById(labelId);
      if (!toggle) return;

      const btns = toggle.querySelectorAll('.format-toggle__btn');
      btns.forEach(btn => {
        btn.addEventListener('click', () => {
          btns.forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          const fmt = btn.dataset.format;
          state[stateKey] = fmt;

          if (fmt === 'position') {
            if (windowGroup) windowGroup.style.display = 'block';
            if (label) {
              label.innerHTML = 'Paste Coordinates <span class="form-sublabel">(Chromosome &lt;tab/space&gt; Position [± Window Applied])</span>';
            }
            if (textarea && !textarea.value.trim()) {
              textarea.placeholder = "14\t1801116\n6\t38024921\n2\t6445100\n20\t31920050\n26\t21150000";
            }
          } else {
            if (windowGroup) windowGroup.style.display = 'none';
            if (label) {
              label.innerHTML = 'Paste Coordinates <span class="form-sublabel">(Chromosome &lt;tab/space&gt; Start &lt;tab/space&gt; End)</span>';
            }
            if (textarea && !textarea.value.trim()) {
              textarea.placeholder = "14\t1790000\t1820000\n6\t37900000\t38100000\n2\t6440000\t6460000\n20\t31880000\t32180000\n26\t21140000\t21170000";
            }
          }
        });
      });
    }

    bindToggle('geneFormatToggle', 'geneWindowGroup', 'geneCoordsInput', 'geneInputLabel', 'geneFormat');
    bindToggle('qtlFormatToggle', 'qtlWindowGroup', 'qtlCoordsInput', 'qtlInputLabel', 'qtlFormat');
    bindToggle('combinedFormatToggle', 'combinedWindowGroup', 'combinedCoordsInput', null, 'combinedFormat');
  }

  // 3. Flanking Window Preset Chips
  function setupWindowChips() {
    function bindChips(containerId, inputId) {
      const container = document.getElementById(containerId);
      const input = document.getElementById(inputId);
      if (!container || !input) return;

      container.querySelectorAll('.chip').forEach(chip => {
        chip.addEventListener('click', () => {
          container.querySelectorAll('.chip').forEach(c => c.classList.remove('active'));
          chip.classList.add('active');
          if (chip.dataset.window) {
            input.value = chip.dataset.window;
          }
        });
      });
    }

    bindChips('geneWindowPresetChips', 'geneCustomWindowInput');
    bindChips('qtlWindowPresetChips', 'qtlCustomWindowInput');
  }

  // 4. Animal Selection Cards Grid
  function setupAnimalGrid() {
    const cards = document.querySelectorAll('#qtlAnimalGrid .animal-select-card');
    const customUploadGroup = document.getElementById('customQtlUploadGroup');

    cards.forEach(card => {
      card.addEventListener('click', () => {
        cards.forEach(c => c.classList.remove('selected'));
        card.classList.add('selected');
        const sp = card.dataset.species;
        state.selectedAnimal = sp;

        if (customUploadGroup) {
          customUploadGroup.style.display = (sp === 'custom') ? 'block' : 'none';
        }

        const qtlInput = document.getElementById('qtlCoordsInput');
        if (qtlInput) {
          if (sp === 'custom') {
            qtlInput.placeholder = "14\t1750000\t1850000\n6\t37900000\t38200000\n2\t6400000\t6500000";
          } else if (ANIMAL_PRESETS[sp]) {
            qtlInput.placeholder = (state.qtlFormat === 'position') ? ANIMAL_PRESETS[sp].pos : ANIMAL_PRESETS[sp].range;
          }
        }

        // Update assembly indicator in results
        const statAssembly = document.getElementById('statAssembly');
        if (statAssembly && card.dataset.assembly) {
          statAssembly.textContent = card.dataset.assembly;
        }
      });
    });
  }

  // 5. Drag & Drop Upload Zones
  function setupUploadZones() {
    function bindZone(zoneId, inputId, infoId, nameId) {
      const zone = document.getElementById(zoneId);
      const input = document.getElementById(inputId);
      const info = document.getElementById(infoId);
      const nameEl = document.getElementById(nameId);
      if (!zone || !input) return;

      zone.addEventListener('click', () => input.click());
      zone.addEventListener('dragover', (e) => { e.preventDefault(); zone.classList.add('drag-over'); });
      zone.addEventListener('dragleave', () => zone.classList.remove('drag-over'));
      zone.addEventListener('drop', (e) => {
        e.preventDefault();
        zone.classList.remove('drag-over');
        if (e.dataTransfer.files.length) handleFile(e.dataTransfer.files[0]);
      });
      input.addEventListener('change', () => {
        if (input.files.length) handleFile(input.files[0]);
      });

      function handleFile(file) {
        if (info && nameEl) {
          nameEl.textContent = `${file.name} (${(file.size / 1024).toFixed(1)} KB)`;
          info.classList.add('visible');
        }
      }
    }

    bindZone('geneCoordsUploadZone', 'geneCoordsFileInput', 'geneCoordsFileInfo', 'geneCoordsFileName');
    bindZone('qtlUploadZone', 'qtlFileInput', 'qtlFileInfo', 'qtlFileName');
    bindZone('combinedUploadZone', 'combinedFileInput', null, null);
  }

  // 6. Example Fill Preset Buttons
  function setupPresetButtons() {
    const geneFillBtn = document.getElementById('loadGeneRangePresetBtn');
    if (geneFillBtn) {
      geneFillBtn.addEventListener('click', () => {
        const input = document.getElementById('geneCoordsInput');
        if (input) {
          input.value = (state.geneFormat === 'position')
            ? "14\t1800000\n6\t38000000\n2\t6445000\n20\t32000000\n26\t21150000\n4\t93270000"
            : "14\t1790000\t1820000\n6\t37900000\t38100000\n2\t6440000\t6460000\n20\t31880000\t32180000\n26\t21140000\t21170000\n4\t93260000\t93290000";
        }
      });
    }

    const qtlFillBtn = document.getElementById('loadQTLRangePresetBtn');
    if (qtlFillBtn) {
      qtlFillBtn.addEventListener('click', () => {
        const input = document.getElementById('qtlCoordsInput');
        const sp = state.selectedAnimal || 'cattle';
        const preset = ANIMAL_PRESETS[sp] || ANIMAL_PRESETS.cattle;
        if (input) {
          input.value = (state.qtlFormat === 'position') ? preset.pos : preset.range;
        }
      });
    }

    const combFillBtn = document.getElementById('loadCombinedPresetBtn');
    if (combFillBtn) {
      combFillBtn.addEventListener('click', () => {
        const input = document.getElementById('combinedCoordsInput');
        if (input) {
          input.value = (state.combinedFormat === 'position')
            ? "14\t1800000\n6\t38000000\n2\t6445000\n20\t32000000"
            : "14\t1750000\t1850000\n6\t37900000\t38200000\n2\t6400000\t6500000\n20\t31800000\t32200000";
        }
      });
    }
  }

  // 7. Execution Handlers
  function setupExecutionButtons() {
    const runGene = document.getElementById('runGeneAnnotation');
    if (runGene) {
      runGene.addEventListener('click', () => {
        const input = document.getElementById('geneCoordsInput');
        if (input && !input.value.trim()) {
          input.value = "14\t1790000\t1820000\n6\t37900000\t38100000\n2\t6440000\t6460000\n20\t31880000\t32180000\n26\t21140000\t21170000";
        }
        executeRun('genes', () => {
          switchResultsTab('genes');
        });
      });
    }

    const runQTL = document.getElementById('runQTLAnnotation');
    if (runQTL) {
      runQTL.addEventListener('click', () => {
        const input = document.getElementById('qtlCoordsInput');
        const sp = state.selectedAnimal || 'cattle';
        const preset = ANIMAL_PRESETS[sp] || ANIMAL_PRESETS.cattle;
        if (input && !input.value.trim()) {
          input.value = preset.range;
        }
        executeRun('qtls', () => {
          switchResultsTab('qtls');
        });
      });
    }

    const runCombined = document.getElementById('runCombinedAnnotation');
    if (runCombined) {
      runCombined.addEventListener('click', () => {
        const input = document.getElementById('combinedCoordsInput');
        if (input && !input.value.trim()) {
          input.value = "14\t1750000\t1850000\n6\t37900000\t38200000\n2\t6400000\t6500000\n20\t31800000\t32200000";
        }
        executeRun('combined', () => {
          switchResultsTab('genes');
        });
      });
    }
  }

  // Example Dataset & Onboarding Logic
  function setupExampleDatasetLoader() {
    const loadExampleBtn = document.getElementById('btnLoadExampleInResults');
    if (loadExampleBtn) {
      loadExampleBtn.addEventListener('click', () => {
        const geneInput = document.getElementById('geneCoordsInput');
        if (geneInput) {
          geneInput.value = "14\t1790000\t1820000\n6\t37900000\t38100000\n2\t6440000\t6460000\n20\t31880000\t32180000\n26\t21140000\t21170000";
        }
        executeRun('genes', () => {
          switchResultsTab('genes');
        });
      });
    }

    // Quick preset buttons under empty state
    document.querySelectorAll('.preset-chip-btn').forEach(btn => {
      btn.addEventListener('click', () => {
        const key = btn.dataset.preset;
        const p = LOCUS_PRESETS[key];
        if (p) {
          // Switch to corresponding mode
          const modeCard = document.getElementById(`type${p.mode.charAt(0).toUpperCase() + p.mode.slice(1)}`);
          if (modeCard) modeCard.click();

          // Set input
          if (p.mode === 'gene') {
            const el = document.getElementById('geneCoordsInput');
            if (el) el.value = p.coords;
          } else if (p.mode === 'qtl') {
            const el = document.getElementById('qtlCoordsInput');
            if (el) el.value = p.coords;
          } else {
            const el = document.getElementById('combinedCoordsInput');
            if (el) el.value = p.coords;
          }

          executeRun(p.mode === 'qtl' ? 'qtls' : 'genes', () => {
            switchResultsTab(p.mode === 'qtl' ? 'qtls' : 'genes');
          });
        }
      });
    });

    const clearBtn = document.getElementById('clearResultsBtn');
    if (clearBtn) {
      clearBtn.addEventListener('click', () => {
        setResultsState(false);
      });
    }
  }

  function setResultsState(hasData) {
    state.hasResults = hasData;
    const emptyBanner = document.getElementById('resultsEmptyBanner');
    const dataContent = document.getElementById('resultsDataContent');

    if (emptyBanner) emptyBanner.style.display = hasData ? 'none' : 'block';
    if (dataContent) dataContent.style.display = hasData ? 'block' : 'none';
  }

  function executeRun(type, callback) {
    const allGenes = (typeof MOCK_GENES !== 'undefined') ? MOCK_GENES : [];
    const allQtls = (typeof MOCK_QTLS !== 'undefined') ? MOCK_QTLS : [];

    state.resultsGenes = allGenes.slice();
    state.resultsQTLs = allQtls.slice();

    renderGenesTable(state.resultsGenes);
    renderQTLsTable(state.resultsQTLs);
    renderTraitChart(state.resultsGenes);

    updateStats();
    setResultsState(true);

    // Scroll smoothly to results
    const resultsEl = document.getElementById('studioResults');
    if (resultsEl) {
      resultsEl.scrollIntoView({ behavior: 'smooth' });
    }

    if (callback) callback();
  }

  function updateStats() {
    const statG = document.getElementById('statGenesCount');
    const statQ = document.getElementById('statQtlsCount');
    const statT = document.getElementById('statTraitsCount');
    const badgeG = document.getElementById('genesCountBadge');
    const badgeQ = document.getElementById('qtlsCountBadge');

    if (statG) statG.textContent = state.resultsGenes.length;
    if (badgeG) badgeG.textContent = state.resultsGenes.length;

    if (statQ) statQ.textContent = state.resultsQTLs.length;
    if (badgeQ) badgeQ.textContent = state.resultsQTLs.length;

    if (statT) {
      const traits = new Set(state.resultsGenes.map(g => g.trait));
      statT.textContent = traits.size;
    }
  }

  // 8. Results Tabs
  function setupResultsTabs() {
    const tabGenes = document.getElementById('resultsTabGenes');
    const tabQTLs = document.getElementById('resultsTabQTLs');
    const tabChart = document.getElementById('resultsTabChart');

    if (tabGenes) tabGenes.addEventListener('click', () => switchResultsTab('genes'));
    if (tabQTLs) tabQTLs.addEventListener('click', () => switchResultsTab('qtls'));
    if (tabChart) tabChart.addEventListener('click', () => switchResultsTab('chart'));
  }

  function switchResultsTab(tab) {
    state.activeResultsTab = tab;
    const tabGenes = document.getElementById('resultsTabGenes');
    const tabQTLs = document.getElementById('resultsTabQTLs');
    const tabChart = document.getElementById('resultsTabChart');

    const viewGenes = document.getElementById('viewGenes');
    const viewQTLs = document.getElementById('viewQTLs');
    const viewChart = document.getElementById('viewChart');

    if (tabGenes) tabGenes.classList.toggle('active', tab === 'genes');
    if (tabQTLs) tabQTLs.classList.toggle('active', tab === 'qtls');
    if (tabChart) tabChart.classList.toggle('active', tab === 'chart');

    if (viewGenes) viewGenes.style.display = (tab === 'genes') ? 'block' : 'none';
    if (viewQTLs) viewQTLs.style.display = (tab === 'qtls') ? 'block' : 'none';
    if (viewChart) {
      viewChart.style.display = (tab === 'chart') ? 'block' : 'none';
      if (tab === 'chart') {
        setTimeout(() => renderTraitChart(state.resultsGenes), 50);
      }
    }
  }

  // Render Genes Table
  function renderGenesTable(genes) {
    const tbody = document.getElementById('genesTableBody');
    const summary = document.getElementById('genesCountSummary');
    if (!tbody) return;

    tbody.innerHTML = '';
    if (summary) summary.textContent = `Showing ${genes.length} genes`;

    if (genes.length === 0) {
      tbody.innerHTML = `<tr><td colspan="9" style="text-align:center; padding: 32px; color: #64748b;">No overlapping genes found.</td></tr>`;
      return;
    }

    genes.forEach((g, idx) => {
      const tr = document.createElement('tr');
      tr.innerHTML = `
        <td>
          <div class="gene-cell">
            <span class="gene-symbol">${g.symbol}</span>
          </div>
        </td>
        <td><span class="gene-id">${g.id}</span></td>
        <td><span class="chr-badge">${g.chr}</span></td>
        <td><span class="coord-tag">${g.start.toLocaleString()} - ${g.end.toLocaleString()}</span></td>
        <td><span style="font-family: var(--font-mono); font-size: 11px; color: #64748b;">${g.biotype || 'protein_coding'}</span></td>
        <td>
          <div class="overlap-meter">
            <div class="overlap-bar" style="width: ${g.overlap_pct}%"></div>
            <span class="overlap-text">${g.overlap_pct}% (${g.overlap_bp.toLocaleString()} bp)</span>
          </div>
        </td>
        <td>
          <span class="tss-tag ${g.distance_to_tss === 0 ? 'tss-zero' : ''}">
            ${g.distance_to_tss === 0 ? 'In-TSS / Exonic' : `+${g.distance_to_tss} bp`}
          </span>
        </td>
        <td><span class="trait-tag">${g.trait}</span></td>
        <td>
          <button class="inspect-btn" data-index="${idx}">
            <svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="11" cy="11" r="8"/><line x1="21" y1="21" x2="16.65" y2="16.65"/></svg>
            Inspect
          </button>
        </td>
      `;
      tbody.appendChild(tr);
    });

    tbody.querySelectorAll('.inspect-btn').forEach(btn => {
      btn.addEventListener('click', () => {
        const i = parseInt(btn.dataset.index, 10);
        openInspectorModal(genes[i]);
      });
    });
  }

  // Render QTLs Table
  function renderQTLsTable(qtls) {
    const tbody = document.getElementById('qtlsTableBody');
    const summary = document.getElementById('qtlsCountSummary');
    if (!tbody) return;

    tbody.innerHTML = '';
    if (summary) summary.textContent = `Showing ${qtls.length} QTLs`;

    if (qtls.length === 0) {
      tbody.innerHTML = `<tr><td colspan="8" style="text-align:center; padding: 32px; color: #64748b;">No overlapping QTLs found.</td></tr>`;
      return;
    }

    qtls.forEach(q => {
      const tr = document.createElement('tr');
      tr.innerHTML = `
        <td><span class="chr-badge" style="background: #faf5ff; color: #7c3aed; border-color: #ddd6fe;">${q.id}</span></td>
        <td><strong style="color: #0f172a;">${q.trait}</strong></td>
        <td><span style="font-family: var(--font-mono); font-size: 11px; color: #b45309;">${q.category}</span></td>
        <td><span class="chr-badge">${q.chr}</span></td>
        <td><span class="coord-tag">${q.start.toLocaleString()} - ${q.end.toLocaleString()}</span></td>
        <td><span style="font-family: var(--font-mono); color: #059669; font-weight: 600;">${q.overlap_pct}%</span></td>
        <td><a href="https://pubmed.ncbi.nlm.nih.gov/${q.pubmed.replace('PMID:', '')}" target="_blank" rel="noopener" style="color: var(--oblo-blue); text-decoration: none; font-family: var(--font-mono); font-size: 11px; font-weight: 600;">${q.pubmed} &nearr;</a></td>
        <td><span style="font-weight: 700; color: #0f172a;">${q.candidate || '-'}</span></td>
      `;
      tbody.appendChild(tr);
    });
  }

  // Render Canvas Trait Bar Chart
  function renderTraitChart(genes) {
    const canvas = document.getElementById('traitChartCanvas');
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    const dpr = window.devicePixelRatio || 1;

    const rect = canvas.getBoundingClientRect();
    canvas.width = rect.width * dpr;
    canvas.height = 240 * dpr;
    ctx.scale(dpr, dpr);

    const w = rect.width;
    const h = 240;

    ctx.clearRect(0, 0, w, h);

    const traitCounts = {};
    genes.forEach(g => {
      traitCounts[g.trait] = (traitCounts[g.trait] || 0) + 1;
    });

    if (Object.keys(traitCounts).length === 0) {
      traitCounts["Milk Fat %"] = 4;
      traitCounts["Fatty Acid Composition"] = 3;
      traitCounts["Milk Yield"] = 2;
      traitCounts["Muscle Development"] = 2;
    }

    const traits = Object.keys(traitCounts);
    const maxVal = Math.max(...Object.values(traitCounts), 1);
    const barWidth = Math.min(64, (w - 80) / traits.length - 20);
    const chartBottom = h - 45;
    const chartHeight = chartBottom - 30;

    // Baseline
    ctx.strokeStyle = '#e2e8f0';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(40, chartBottom);
    ctx.lineTo(w - 20, chartBottom);
    ctx.stroke();

    const colors = ['#f05a42', '#f4c03b', '#113cbd', '#db7ea6', '#3b82f6', '#10b981'];

    traits.forEach((t, i) => {
      const count = traitCounts[t];
      const barH = (count / maxVal) * chartHeight;
      const x = 50 + i * (barWidth + 28);
      const y = chartBottom - barH;

      const grad = ctx.createLinearGradient(x, y, x, chartBottom);
      grad.addColorStop(0, colors[i % colors.length]);
      grad.addColorStop(1, 'rgba(25, 70, 230, 0.08)');

      ctx.fillStyle = grad;
      ctx.beginPath();
      ctx.roundRect ? ctx.roundRect(x, y, barWidth, barH, [6, 6, 0, 0]) : ctx.rect(x, y, barWidth, barH);
      ctx.fill();

      // Count label
      ctx.fillStyle = '#0f172a';
      ctx.font = '700 12px JetBrains Mono, monospace';
      ctx.textAlign = 'center';
      ctx.fillText(count.toString(), x + barWidth / 2, y - 8);

      // Trait label
      ctx.fillStyle = '#475569';
      ctx.font = '500 11px Plus Jakarta Sans, sans-serif';
      ctx.textAlign = 'center';
      const shortName = t.length > 16 ? t.slice(0, 14) + '…' : t;
      ctx.fillText(shortName, x + barWidth / 2, chartBottom + 22);
    });
  }

  // 9. Search Filters
  function setupSearchFilters() {
    const genesSearch = document.getElementById('genesTableSearch');
    if (genesSearch) {
      genesSearch.addEventListener('input', (e) => {
        const q = e.target.value.toLowerCase();
        const filtered = state.resultsGenes.filter(g =>
          g.symbol.toLowerCase().includes(q) ||
          g.id.toLowerCase().includes(q) ||
          g.chr.toLowerCase().includes(q) ||
          g.trait.toLowerCase().includes(q)
        );
        renderGenesTable(filtered);
      });
    }

    const qtlsSearch = document.getElementById('qtlsTableSearch');
    if (qtlsSearch) {
      qtlsSearch.addEventListener('input', (e) => {
        const q = e.target.value.toLowerCase();
        const filtered = state.resultsQTLs.filter(item =>
          item.id.toLowerCase().includes(q) ||
          item.trait.toLowerCase().includes(q) ||
          item.category.toLowerCase().includes(q) ||
          item.chr.toLowerCase().includes(q)
        );
        renderQTLsTable(filtered);
      });
    }
  }

  // 10. Comprehensive Download Suite (TSV, CSV, JSON, Report)
  function setupDownloadButtons() {
    const btnTsv = document.getElementById('downloadTsvBtn');
    const btnCsv = document.getElementById('downloadCsvBtn');
    const btnJson = document.getElementById('downloadJsonBtn');
    const btnReport = document.getElementById('downloadReportBtn');

    if (btnTsv) {
      btnTsv.addEventListener('click', () => exportDelimited('\t', 'qgat_results.tsv'));
    }
    if (btnCsv) {
      btnCsv.addEventListener('click', () => exportDelimited(',', 'qgat_results.csv'));
    }
    if (btnJson) {
      btnJson.addEventListener('click', () => exportJson('qgat_results.json'));
    }
    if (btnReport) {
      btnReport.addEventListener('click', () => exportReport('qgat_summary_report.txt'));
    }
  }

  function exportDelimited(delimiter, filename) {
    const isGenes = state.activeResultsTab === 'genes' || state.activeResultsTab === 'chart';
    let content = '';

    if (isGenes) {
      const headers = ['Symbol', 'Gene_ID', 'Chromosome', 'Start', 'End', 'Biotype', 'Overlap_bp', 'Overlap_pct', 'TSS_Distance', 'Trait_Association', 'Key_Variant'];
      const rows = state.resultsGenes.map(g => [
        g.symbol,
        g.id,
        g.chr,
        g.start,
        g.end,
        g.biotype,
        g.overlap_bp,
        g.overlap_pct,
        g.distance_to_tss,
        `"${g.trait}"`,
        `"${g.key_variant || ''}"`
      ]);
      content = [headers.join(delimiter), ...rows.map(r => r.join(delimiter))].join('\n');
    } else {
      const headers = ['QTL_ID', 'Trait', 'Category', 'Chromosome', 'Start', 'End', 'Overlap_pct', 'PubMed', 'Candidate_Gene'];
      const rows = state.resultsQTLs.map(q => [
        q.id,
        `"${q.trait}"`,
        `"${q.category}"`,
        q.chr,
        q.start,
        q.end,
        q.overlap_pct,
        q.pubmed,
        q.candidate || ''
      ]);
      content = [headers.join(delimiter), ...rows.map(r => r.join(delimiter))].join('\n');
    }

    triggerDownload(content, filename, 'text/plain;charset=utf-8;');
  }

  function exportJson(filename) {
    const payload = {
      generated_by: "QGAT (Quantitative Trait Loci & Gene Annotation Tool)",
      date: new Date().toISOString(),
      species: state.selectedAnimal,
      genes_count: state.resultsGenes.length,
      qtls_count: state.resultsQTLs.length,
      genes: state.resultsGenes,
      qtls: state.resultsQTLs
    };
    triggerDownload(JSON.stringify(payload, null, 2), filename, 'application/json');
  }

  function exportReport(filename) {
    const dateStr = new Date().toLocaleString();
    const text = [
      "============================================================",
      "             QGAT ANNOTATION SUMMARY REPORT                 ",
      "        Quantitative Trait Loci & Gene Annotation Tool      ",
      "============================================================",
      `Generated: ${dateStr}`,
      `Target Species: ${state.selectedAnimal}`,
      `Total Genes Annotated: ${state.resultsGenes.length}`,
      `Total Overlapping QTLs: ${state.resultsQTLs.length}`,
      "------------------------------------------------------------",
      "TOP CANDIDATE GENES DISCOVERED:",
      ...state.resultsGenes.map((g, i) =>
        `${i + 1}. ${g.symbol} (${g.id}) | ${g.chr}:${g.start}-${g.end} | Overlap: ${g.overlap_pct}% | Trait: ${g.trait}\n   Causal Variant: ${g.key_variant || 'None documented'}`
      ),
      "------------------------------------------------------------",
      "ANIMAL QTLdb TRAIT INTERVALS:",
      ...state.resultsQTLs.map((q, i) =>
        `${i + 1}. [${q.id}] ${q.trait} (${q.category}) | ${q.chr}:${q.start}-${q.end} | Overlap: ${q.overlap_pct}% | ${q.pubmed}`
      ),
      "============================================================",
      "For details, visit: https://akanksha-ch0.github.io/QGAT/"
    ].join('\n');

    triggerDownload(text, filename, 'text/plain;charset=utf-8;');
  }

  function triggerDownload(content, filename, mimeType) {
    const blob = new Blob([content], { type: mimeType });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }

  // 11. Modal / Exon Inspector
  function setupModal() {
    const modal = document.getElementById('inspectorModal');
    const closeBtn = document.getElementById('modalCloseBtn');
    if (!modal) return;

    if (closeBtn) {
      closeBtn.addEventListener('click', () => modal.classList.remove('active'));
    }
    modal.addEventListener('click', (e) => {
      if (e.target === modal) modal.classList.remove('active');
    });
  }

  function openInspectorModal(gene) {
    const modal = document.getElementById('inspectorModal');
    if (!modal || !gene) return;

    document.getElementById('modalGeneSymbol').textContent = gene.symbol;
    document.getElementById('modalGeneId').textContent = gene.id;
    document.getElementById('modalTrait').textContent = gene.trait || 'Candidate Locus';
    document.getElementById('modalChr').textContent = `${gene.chr}:${gene.start.toLocaleString()}-${gene.end.toLocaleString()}`;
    document.getElementById('modalFunction').textContent = gene.function || 'Functional annotation pending.';
    document.getElementById('modalVariant').textContent = gene.key_variant || 'No causal mutation recorded.';

    const exonTrack = document.getElementById('modalExonTrack');
    if (exonTrack && gene.exons) {
      exonTrack.innerHTML = '';
      gene.exons.forEach(ex => {
        const box = document.createElement('div');
        box.className = `exon-block ${ex.is_causal ? 'causal' : ''} ${ex.type}`;
        box.title = `Exon ${ex.num}: ${ex.len} bp (${ex.type}) ${ex.variant_note ? '• ' + ex.variant_note : ''}`;
        box.style.flex = ex.len;
        box.innerHTML = `<span class="exon-num">${ex.num}</span>`;
        exonTrack.appendChild(box);
      });
    }

    const pathwayContainer = document.getElementById('modalPathways');
    if (pathwayContainer) {
      pathwayContainer.innerHTML = '';
      if (gene.pathways && gene.pathways.length > 0) {
        gene.pathways.forEach(p => {
          const tag = document.createElement('span');
          tag.className = 'pathway-tag';
          tag.textContent = p;
          pathwayContainer.appendChild(tag);
        });
      }
    }

    modal.classList.add('active');
  }

  // 12. CLI Tabs & Copy
  function setupCliTabs() {
    const tabBtns = document.querySelectorAll('.cli-tab-btn');
    const tabPanes = document.querySelectorAll('.cli-code-block');

    tabBtns.forEach(btn => {
      btn.addEventListener('click', () => {
        tabBtns.forEach(b => b.classList.remove('active'));
        tabPanes.forEach(p => p.classList.remove('active'));

        btn.classList.add('active');
        const target = btn.dataset.target;
        const pane = document.getElementById(target);
        if (pane) pane.classList.add('active');
      });
    });

    const copyBtns = document.querySelectorAll('.copy-cli-btn');
    copyBtns.forEach(btn => {
      btn.addEventListener('click', () => {
        const text = btn.dataset.copy || '';
        navigator.clipboard.writeText(text).then(() => {
          const original = btn.innerHTML;
          btn.innerHTML = `<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5"><polyline points="20 6 9 17 4 12"/></svg> Copied!`;
          btn.classList.add('copied');
          setTimeout(() => {
            btn.innerHTML = original;
            btn.classList.remove('copied');
          }, 2000);
        });
      });
    });
  }

  // Boot on DOM ready
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }

})();
