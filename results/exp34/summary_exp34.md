# Exp34 — Ultra iyon itki denklemi: bütün veriden (Exp0–Exp33) optimize edilmiş çalışma zarfı

Tarih: 2026-10-16. Kod: `experiments/exp34_ultra_thrust/ultra34.py`; şekil kodu `fig34.py`.
Öngörüler optimizasyon ızgarası hesaplanmadan önce yazıldı: `predictions_exp34.json` (UT1–UT6; dayanak kör-olmayan olarak bildirildi).
Sonuç: `thrust34.json` (optima, Pareto cepheleri, hücre ızgaraları, görevler, doğrulama). Şekil: `fig_exp34.png` (4500×3000 px) — Pareto cephesi, tür optimumları, hücre ızgarası ve referans görev panelleri.

## Denklem (nihai biçim)

```
F        = gamma_b * I_tot * sqrt(2 m V / q),        I_tot = N * I_c,   I_c = 0.75 * j_on * w_s(lambda, r_D)
P        = I_tot * (V + V_M + E_giris/e) / eta_rec
F/P      = 2 eta / v_e,    eta = gamma_b * V / (V + V_M + E_giris/e),   v_e = sqrt(2 q V / m),   Isp = v_e / g0
Delta_v  = v_e ln(M0/M1)   [yerleşik yakıt]   |   v_e = v_orb (1 + C_D/(2 eps f_m))  [ortam atomu, Exp28 EX3]
F/P_isin = eta_pv * F/P    [ışınlanan güç, EX4]
```

Kısıtlar (hepsi proje verisinden):
- uzay yükü: `F/A <= (8/9) eps0 (V/d)^2`, `J <= J_CL(V,d)`;
- yayma beslemesi: `beta * V/d >= F_th(tür)` (Exp24 bariyeri: önceden yüklü 1,97 V/nm; In buharlaşma 8,26; Cs/Xe WKB 3,80/20,89) ve hız tarafı `I_tot/q <= k(F,T) * N_aktif_bölge`, `k = nu exp(-(Q0 - alpha F)/kT)` (Exp24) — hız tarafı bu deneyde sayısal kapatılmadı, açık maddedir;
- geometri tabanı `d >= 4,148 nm` (SKY130 t_ox);
- hücre basamakları: N tam sayı (Exp19/20: k = 1, 2, 4, zorla 3→7), I_c gövde-etkisi bandında (aşağıda);
- foton kapısı yalnız yelkende 2/c (Exp27); kapalı araçta F/P ≤ 2η/v_e her zaman.

gamma_b = 0,9 (Exp28/32 demet katsayısı), eta_rec nominal 1,0 (Exp23 aktarım aralığı 0,985–1,023). Görev yanma süreleri hem model eta hem Exp28'in eta = 0,6 referansıyla verilir (PPU/ızgara kayıpları model eta'da yoktur).

## Hükümler (UT1–UT6, kayıtlı ölçütlere göre)

- **UT1 ✓ (iç optimum).** F/P(V) her türde V* = C = V_M + E_giris/e'de tepe yapar (ızgara = analitik, <%0,1). EMI⁺: V* = 1,70 V, **F/P = 523,8 µN/W**, Isp = 175 s. In⁺: 272,6 µN/W @6,49 V. Tür sıralaması (kısıtsız tepe): EMI⁺ 524 > Cs⁺ 348 > In⁺ 273 > Ga⁺ 209 > Xe⁺ 207 > O⁻ 199 µN/W.
- **UT2 kısmen.** EMI⁺ fizibil optimum yayma sınırında: V = 8,36 V (sınır 8,18), **F/P = 392,7 µN/W**, Isp = 388 s ✓ kayıtlı aralıkta. In sınırı 34,25 V, orada F/P = 198,1 µN/W ✓. **Falsifiye olan alt madde:** "uç beta ≥ 50 gerekir" — hesaplanan beta (V*'ta, d_min'de) yalnızca **5,28** (EMI⁺'da 4,81). Yani makul bir uç yükseltimi (beta ≈ 5) kısıtsız optimumu bile fizibil kılar; bu denklem için iyi haberdir ama kayıtlı sayı yanlıştı.
- **UT3 ✗ (kayıtlı aralık yanlış; fizik sonucu daha güçlü).** Yayma sınırındaki CL itki yoğunluğu kayıtlı [2,5–3,5]·10⁴ değil, **3,06·10⁷ N/m²** çıktı (kayıtta 10³ kat aritmetik kayması). Buna bağlı P/A = 7,80·10¹⁰ W/m² ve 100 W'da CL-doygun alan 1,28·10⁻⁹ m² (~36×36 µm). Alt iddialar (NSTAR'ın ≥10⁴ katı — gerçekte 2,4·10⁷ kat; tasarım güç/ısıl sınırlıdır, yoğunluk sınırlı değildir) daha güçlü biçimde tutar, fakat hüküm kayıtlı ölçüte göre **fail**'dir.
- **UT4 kısmen.** Δv tarafı tuttu: In @5 kV Δv = 20,46 km/s (yanma eta=0,6 ile 4,44 yıl), nano-aralık optimumu aynı yakıtla yalnız 0,85 km/s. **Falsifiye olan alt madde:** diz (Δv ≥ 10 km/s koşuluyla en yüksek F/P) Isp'si kayıtlı [2000, 4500] s aralığında değil, tanım gereği **4570 s**'dir (v_e = 10/ln(10/8) = 44,81 km/s). Diz noktası: In @1195 V, **F/P = 39,9 µN/W** (≈ NSTAR'ın ölçülen 40 µN/W'si, Isp %47 daha yüksekte).
- **UT5 ✗ (kayıtlı biçimiyle; çekirdek madde tuttu).** s = 0,33 sabitleriyle 9 ızgara noktasından tam biri bölünür (r_D20/λ120: λ_c = 106 µm, w_s = 519 µm, I_c = 0,33 mA ✓ kayıtlı aralıkta). Ama "s = 0,50'de dokuzu da bölünür" tutmadı: **6/9** (r_D5'te λ40/60 ve r_D10'da λ40 bölünmez; λ_c = 68,7/48,6 µm). Hücre beslemesi bandı (bölünen noktalarda): s=1: 0,99–2,60 mA; s=0,50: 0,40–1,00 mA; s=0,33: 0,33 mA (tek köşe). 5 kV In'da hücre başına (λ60/r_D10): 149 µN (s=1), 96 µN (s=0,50), yok (s=0,33).
- **UT6 ✓ (gerileme yok).** Çekirdek Exp29 ile cebirsel olarak aynı olduğundan NEXT kapanışı aynen korunur (kayıtlı MC1: 40/40, ortalama +%1,01, RMS %1,45); beş ölçülmüş itkiç satırında eta yeniden hesabı saklı değerlerle birebir (sapma 0). Ayrıca F/P(In, 5 kV) = 19,61 µN/W, aer28 `cells` çıktısıyla aynı.

## Çalışma zarfı (rejimler)

- **R1 — nano-aralık, önceden yüklü (maks F/P):** EMI⁺ @8,4 V (beta=1): 392,7 µN/W, Isp 388 s; beta ≥ 4,8 uçla @1,7 V: 523,8 µN/W, Isp 175 s. Foton sınırının (3,34 nN/W) 1,2–1,6·10⁵ katı, NSTAR F/P'sinin ~10 katı — bedeli düşük Isp'dir (Δv ≤ 0,85 km/s, 2 kg). Kullanım: kısa süreli yüksek itki, istasyon tutma, sürüklenme telafisi. Isıl not: CL yoğunluğuna çıkılmaz; alan, güç yoğunluğu taşınabilir düzeyde seçilir, F/P denklemi değişmez.
- **R2 — hücreli kV (görev):** In⁺ @1,2–5 kV: F/P 39,9 → 19,6 µN/W, Isp 4570 → 9348 s, Δv 10 → 20,5 km/s. İtki, I_c basamaklarıyla (k = 1, 2, 4, 7…) kuantumlu ayarlanır; hücre başına 149 µN @5 kV (s=1).
- **R3 — ortam soluma (200 km):** bütün türler toplanırsa v_e = 29,21 km/s ve 443 W/m² ağız ile sürüklenme karşılanır, Δv yerleşik kütle oranıyla sınırlanmaz; yalnız O ile 72,3 km/s ve 901 W/m² (toplam yoğunluk 2·10⁻¹⁰ kg/m³ varsayımı).
- **R4 — ışınlanan güç + atom:** eta_pv = 0,5 ile nano noktada ışınlanan watt başına 196 µN, yelkenin (2/c) **29.431 katı**; diz noktasında 2.994 katı (EX4 ile uyumlu).
- Pareto (EMI⁺, beta=1, s=1): Isp 388→500→1000→2000→3000→4570→9000 s için F/P 393→327→178→90→60→40→20 µN/W. Gövde bandı (s=0,50/0,33) R1 tepesini yalnız 524→491→466 µN/W oynatır; belirsizlik esas olarak R2'nin I_c beslemesindedir.

## Açık maddeler (denklemi kapatmayanlar)

1. **Yayma hızı tarafı:** k(F,T)·N_bölge besleme akımı, CL talebiyle sayısal karşılaştırılmadı; nano-aralık noktasında gerekli bölge yoğunluğu ve ısınma Exp35 adayıdır.
2. **İleri gövde etkisi** (Exp33 C33-1/2): I_c bandı 0,33–2,60 mA; laboratuvarın ilk işi ileri-VB ölçümüdür.
3. **Fiziksel doğrulama yok:** R1'in EMI⁺/bariyer sayıları Exp24 model yasası + literatür çapalarıdır; R2'nin hücre kilidi/histerezisi Exp33 protokolüyle laboratuvar bekler. Denklemin itki çekirdeği (γ√·) Exp29/30 ile ölçüm-kapalıdır; hücre ve nano-aralık terimleri ölçüme açıktır.

## Kaynaklar

Exp28 denklemi §6 ve `aer28.py`/`explore28.py` (yeniden koşuda birebir doğrulandı), Exp29 `closure.json`, Exp32 `nist_sky130.json` sabitleri, Exp33 `cell33_sensitivity.json` (s-bantları, keşif), Exp24 bariyer yasası, Exp17 `law.py`, digest_exp0_32.md.
