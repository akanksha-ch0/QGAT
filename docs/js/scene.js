/**
 * QGAT — OBLO 3D Interactive Clay Scene
 * Crafted with Three.js inspired by ThreeUI Oblo Landing Page aesthetic
 */

(function () {
  'use strict';

  const container = document.getElementById('obloCanvasContainer');
  const canvas = document.getElementById('obloCanvas');
  if (!container || !canvas || typeof THREE === 'undefined') {
    return;
  }

  // Scene & Camera
  const scene = new THREE.Scene();

  const camera = new THREE.PerspectiveCamera(
    38,
    container.clientWidth / container.clientHeight,
    0.1,
    100
  );
  camera.position.set(0, 1.1, 8.8);
  camera.lookAt(0, 0.6, 0);

  // Renderer
  const renderer = new THREE.WebGLRenderer({
    canvas: canvas,
    alpha: true,
    antialias: true,
    powerPreference: 'high-performance'
  });
  renderer.setSize(container.clientWidth, container.clientHeight);
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  renderer.shadowMap.enabled = true;
  renderer.shadowMap.type = THREE.PCFSoftShadowMap;
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 1.05;

  // Master Group for mouse parallax
  const masterGroup = new THREE.Group();
  scene.add(masterGroup);

  // Lighting
  const ambientLight = new THREE.HemisphereLight(0xffffff, 0x142870, 0.85);
  scene.add(ambientLight);

  const sunLight = new THREE.DirectionalLight(0xfffaed, 1.4);
  sunLight.position.set(6, 9, 6);
  sunLight.castShadow = true;
  sunLight.shadow.mapSize.width = 1024;
  sunLight.shadow.mapSize.height = 1024;
  sunLight.shadow.camera.near = 1;
  sunLight.shadow.camera.far = 20;
  sunLight.shadow.bias = -0.0005;
  scene.add(sunLight);

  const fillLight = new THREE.DirectionalLight(0x8cb6ff, 0.7);
  fillLight.position.set(-6, 4, 3);
  scene.add(fillLight);

  const rimLight = new THREE.DirectionalLight(0xffeedd, 0.5);
  rimLight.position.set(0, -3, -5);
  scene.add(rimLight);

  // Clay Material Creator
  function createClay(colorHex, roughness = 0.86, metalness = 0.03) {
    return new THREE.MeshStandardMaterial({
      color: colorHex,
      roughness: roughness,
      metalness: metalness,
      flatShading: false
    });
  }

  // Color Palette (Signature Oblo Clay)
  const CLAY_CORAL = 0xf05a42;      // Arch & accent ball
  const CLAY_COBALT = 0x113cbd;     // Stairs & base sphere
  const CLAY_WHITE = 0xf3f2ec;      // Pedestal & Torus ring
  const CLAY_OCHRE = 0xf4c03b;      // Cone
  const CLAY_PINK = 0xdb7ea6;       // Quarter cylinder & pyramid
  const CLAY_YELLOW_DISC = 0xf4cc4c;// Upright semicircular disc
  const CLAY_SKY = 0x5a88f5;        // Floating high ball

  const interactives = [];

  function registerInteractive(mesh, name, bounceAmount = 0.25) {
    mesh.castShadow = true;
    mesh.receiveShadow = true;
    mesh.userData = {
      name: name,
      baseY: mesh.position.y,
      baseRotX: mesh.rotation.x,
      baseRotY: mesh.rotation.y,
      baseRotZ: mesh.rotation.z,
      bounce: 0,
      bounceVel: 0,
      hoverScale: 1
    };
    interactives.push(mesh);
    return mesh;
  }

  // 1. Pedestal (Chalk White Platform)
  const pedestalGeo = new THREE.CylinderGeometry(3.6, 3.8, 0.45, 64);
  const pedestalMat = createClay(CLAY_WHITE, 0.92, 0.02);
  const pedestal = new THREE.Mesh(pedestalGeo, pedestalMat);
  pedestal.position.set(0, -0.22, 0);
  pedestal.receiveShadow = true;
  masterGroup.add(pedestal);

  // Pedestal rim ring
  const rimGeo = new THREE.TorusGeometry(3.68, 0.06, 16, 64);
  rimGeo.rotateX(Math.PI / 2);
  const rim = new THREE.Mesh(rimGeo, pedestalMat);
  rim.position.set(0, -0.02, 0);
  masterGroup.add(rim);

  // 2. The Coral Arch with center aperture
  (function createCoralArch() {
    const shape = new THREE.Shape();
    const w = 2.4;
    const h = 3.2;
    const r = w / 2;

    // Outer contour with rounded arch top
    shape.moveTo(-r, 0);
    shape.lineTo(-r, h - r);
    shape.absarc(0, h - r, r, Math.PI, 0, true);
    shape.lineTo(r, 0);
    shape.lineTo(-r, 0);

    // Inner aperture hole
    const hole = new THREE.Path();
    const hw = 1.1;
    const hh = 2.1;
    const hr = hw / 2;
    hole.moveTo(-hr, 0);
    hole.lineTo(-hr, hh - hr);
    hole.absarc(0, hh - hr, hr, Math.PI, 0, true);
    hole.lineTo(hr, 0);
    hole.lineTo(-hr, 0);
    shape.holes.push(hole);

    const extrudeSettings = {
      depth: 0.85,
      bevelEnabled: true,
      bevelSegments: 4,
      steps: 1,
      bevelSize: 0.05,
      bevelThickness: 0.05
    };

    const archGeo = new THREE.ExtrudeGeometry(shape, extrudeSettings);
    archGeo.center();

    const archMat = createClay(CLAY_CORAL, 0.85);
    const arch = new THREE.Mesh(archGeo, archMat);
    arch.position.set(0.1, 1.65, -0.15);
    registerInteractive(arch, 'Coral Arch');
    masterGroup.add(arch);

    // Dome on top of arch (Pearl white half-dome)
    const domeGeo = new THREE.SphereGeometry(0.55, 32, 16, 0, Math.PI * 2, 0, Math.PI / 2);
    const domeMat = createClay(CLAY_WHITE, 0.9);
    const dome = new THREE.Mesh(domeGeo, domeMat);
    dome.position.set(0.1, 3.25, 0.28);
    registerInteractive(dome, 'Summit Dome');
    masterGroup.add(dome);
  })();

  // 3. Cobalt Geometric Stairs (Stepped Pyramid Block)
  (function createCobaltStairs() {
    const stairsGroup = new THREE.Group();
    stairsGroup.position.set(1.45, 0.0, 0.1);

    const stepMat = createClay(CLAY_COBALT, 0.84);
    const numSteps = 4;
    const totalW = 1.25;
    const totalD = 1.1;
    const totalH = 1.75;
    const stepH = totalH / numSteps;

    for (let i = 0; i < numSteps; i++) {
      const stepD = totalD * (1 - (i * 0.22));
      const boxGeo = new THREE.BoxGeometry(totalW, stepH, stepD);
      const box = new THREE.Mesh(boxGeo, stepMat);
      box.position.set(0, stepH * (i + 0.5), (totalD - stepD) / 2);
      box.castShadow = true;
      box.receiveShadow = true;
      stairsGroup.add(box);
    }

    registerInteractive(stairsGroup, 'Cobalt Steps');
    masterGroup.add(stairsGroup);
  })();

  // 4. White Chalk Torus Ring (leaning forward)
  (function createWhiteTorus() {
    const torusGeo = new THREE.TorusGeometry(0.95, 0.32, 32, 64);
    const torusMat = createClay(CLAY_WHITE, 0.88);
    const torus = new THREE.Mesh(torusGeo, torusMat);
    torus.position.set(1.7, 0.85, 0.95);
    torus.rotation.set(0.25, -0.4, 0.45);
    registerInteractive(torus, 'Chalk Torus');
    masterGroup.add(torus);
  })();

  // 5. Warm Ochre Cone
  (function createOchreCone() {
    const coneGeo = new THREE.ConeGeometry(0.55, 1.35, 48);
    const coneMat = createClay(CLAY_OCHRE, 0.85);
    const cone = new THREE.Mesh(coneGeo, coneMat);
    cone.position.set(0.7, 0.68, 1.25);
    registerInteractive(cone, 'Ochre Cone');
    masterGroup.add(cone);
  })();

  // 6. Pink Quarter Cylinder Block
  (function createPinkQuarter() {
    const quarterGeo = new THREE.CylinderGeometry(1.2, 1.2, 0.9, 32, 1, false, 0, Math.PI / 2);
    const quarterMat = createClay(CLAY_PINK, 0.88);
    const quarter = new THREE.Mesh(quarterGeo, quarterMat);
    quarter.position.set(-1.1, 0.45, 0.35);
    quarter.rotation.set(0, -Math.PI / 3, 0);
    registerInteractive(quarter, 'Curved Quarter');
    masterGroup.add(quarter);
  })();

  // 7. Yellow Semicircular Backing Disc (Upright disc behind arch)
  (function createYellowDisc() {
    const discGeo = new THREE.CylinderGeometry(1.35, 1.35, 0.18, 48, 1, false, 0, Math.PI);
    const discMat = createClay(CLAY_YELLOW_DISC, 0.88);
    const disc = new THREE.Mesh(discGeo, discMat);
    disc.position.set(1.0, 1.75, -0.75);
    disc.rotation.set(0, 0, -Math.PI / 2);
    registerInteractive(disc, 'Sun Disc');
    masterGroup.add(disc);
  })();

  // 8. Floating Spheres & Geometries (with gentle sinusoidal bobbing)
  const floaters = [];

  function addFloater(geo, mat, x, y, z, name, speed = 1.0, amp = 0.12, rotSpeed = 0.5) {
    const mesh = new THREE.Mesh(geo, mat);
    mesh.position.set(x, y, z);
    registerInteractive(mesh, name);
    mesh.userData.floatSpeed = speed;
    mesh.userData.floatAmp = amp;
    mesh.userData.rotSpeed = rotSpeed;
    mesh.userData.phase = Math.random() * Math.PI * 2;
    masterGroup.add(mesh);
    floaters.push(mesh);
    return mesh;
  }

  // Floating Coral Ball (top left)
  addFloater(
    new THREE.SphereGeometry(0.32, 32, 32),
    createClay(CLAY_CORAL, 0.86),
    -1.7, 2.6, 0.6,
    'Terracotta Orb',
    1.2, 0.15, 0.8
  );

  // Small Dark Cobalt Ball (resting on pedestal front)
  addFloater(
    new THREE.SphereGeometry(0.22, 32, 32),
    createClay(CLAY_COBALT, 0.85),
    -0.35, 0.25, 1.6,
    'Cobalt Sphere',
    0.6, 0.04, 0.3
  );

  // Floating Pearl White Ball (right)
  addFloater(
    new THREE.SphereGeometry(0.24, 32, 32),
    createClay(CLAY_WHITE, 0.88),
    2.6, 1.9, 0.4,
    'Pearl Orb',
    1.4, 0.12, 1.0
  );

  // Floating Pale Blue Ball (high top right)
  addFloater(
    new THREE.SphereGeometry(0.18, 32, 32),
    createClay(CLAY_SKY, 0.86),
    2.2, 3.5, -0.6,
    'Sky Orb',
    1.0, 0.14, 0.6
  );

  // Floating Pink Pyramid (top right)
  addFloater(
    new THREE.TetrahedronGeometry(0.36),
    createClay(CLAY_PINK, 0.85),
    2.5, 2.7, 0.1,
    'Pink Pyramid',
    1.6, 0.16, 1.2
  );

  // Raycasting & Pointer Tracking
  const raycaster = new THREE.Raycaster();
  const mouse = new THREE.Vector2(-999, -999);
  let targetRotY = 0;
  let targetRotX = 0;
  let isHoveringCanvas = false;

  function updatePointer(e) {
    const rect = canvas.getBoundingClientRect();
    mouse.x = ((e.clientX - rect.left) / rect.width) * 2 - 1;
    mouse.y = -((e.clientY - rect.top) / rect.height) * 2 + 1;

    // Gentle parallax calculation
    targetRotY = (mouse.x * 0.28);
    targetRotX = -(mouse.y * 0.16);
  }

  window.addEventListener('mousemove', (e) => {
    updatePointer(e);
  });

  // Touch support
  window.addEventListener('touchmove', (e) => {
    if (e.touches.length > 0) {
      updatePointer(e.touches[0]);
    }
  }, { passive: true });

  // Click bounce reaction
  function triggerBounce(mesh) {
    mesh.userData.bounceVel = 0.18;
  }

  canvas.addEventListener('click', (e) => {
    const rect = canvas.getBoundingClientRect();
    const clickMouse = new THREE.Vector2(
      ((e.clientX - rect.left) / rect.width) * 2 - 1,
      -((e.clientY - rect.top) / rect.height) * 2 + 1
    );
    raycaster.setFromCamera(clickMouse, camera);
    const intersects = raycaster.intersectObjects(interactives, true);

    if (intersects.length > 0) {
      // Find top interactive parent
      let obj = intersects[0].object;
      while (obj && !obj.userData.name && obj.parent) {
        obj = obj.parent;
      }
      if (obj && obj.userData) {
        triggerBounce(obj);
      }
    }
  });

  // Animation Loop
  let clock = new THREE.Clock();

  function animate() {
    requestAnimationFrame(animate);

    const elapsedTime = clock.getElapsedTime();

    // Smooth Parallax Damping
    masterGroup.rotation.y += (targetRotY - masterGroup.rotation.y) * 0.05;
    masterGroup.rotation.x += (targetRotX - masterGroup.rotation.x) * 0.05;

    // Floating animation
    for (let i = 0; i < floaters.length; i++) {
      const f = floaters[i];
      const speed = f.userData.floatSpeed || 1.0;
      const amp = f.userData.floatAmp || 0.1;
      const phase = f.userData.phase || 0;
      const rot = f.userData.rotSpeed || 0.5;

      f.position.y = f.userData.baseY + Math.sin(elapsedTime * speed + phase) * amp;
      f.rotation.y = f.userData.baseRotY + elapsedTime * rot * 0.4;
      f.rotation.x = f.userData.baseRotX + Math.sin(elapsedTime * speed * 0.5) * 0.1;
    }

    // Spring physics on interactive bounces
    for (let i = 0; i < interactives.length; i++) {
      const item = interactives[i];
      if (item.userData.bounceVel !== 0 || item.userData.bounce !== 0) {
        item.userData.bounce += item.userData.bounceVel;
        item.userData.bounceVel -= item.userData.bounce * 0.18; // spring constant
        item.userData.bounceVel *= 0.85; // damping

        item.position.y = item.userData.baseY + item.userData.bounce;

        if (Math.abs(item.userData.bounce) < 0.001 && Math.abs(item.userData.bounceVel) < 0.001) {
          item.userData.bounce = 0;
          item.userData.bounceVel = 0;
          item.position.y = item.userData.baseY;
        }
      }
    }

    renderer.render(scene, camera);
  }

  animate();

  // Resize Handler
  function handleResize() {
    if (!container) return;
    const w = container.clientWidth;
    const h = container.clientHeight;
    camera.aspect = w / h;

    // Responsive camera zoom adjustment
    if (w < 768) {
      camera.fov = 48;
      camera.position.z = 11.2;
      camera.position.y = 1.3;
    } else if (w < 1100) {
      camera.fov = 42;
      camera.position.z = 9.8;
      camera.position.y = 1.2;
    } else {
      camera.fov = 38;
      camera.position.z = 8.8;
      camera.position.y = 1.1;
    }
    camera.updateProjectionMatrix();
    renderer.setSize(w, h);
  }

  window.addEventListener('resize', handleResize);
  handleResize();

})();
