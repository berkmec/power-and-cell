# Exp29 — Atom + pozitif enerji denkleminin ölçülmüş veriyle kapanış testi

Tarih: 2026-10-16. Kod: `experiments/exp29_measured_closure/mc29.py` (`python mc29.py fig`). Öngörüler hesaptan önce yazıldı:
`predictions_exp29.json` (14:44). Sonuçlar: `closure.json`. Veri: `next_diamant2010_table2.json` (Diamant 2010 Tablo 2, PDF
metninden programla ayrıştırıldı, 40/40 satır). Şekil: `fig_exp29.png`.

**Neden.** Exp28 denklemi model düzeyindeydi; Stop kancası ölçümle doğrulama istedi. Laboratuvar burada yok. Ama denklemin her
halkası, girdisi ve çıktısı ayrı ayrı ölçülmüş yayımlanmış deneylerle serbest parametre olmadan sınanabilir. Bu deney bunu yapıyor.

**Açıklama (kör değil).** Kaynakları okurken üç NEXT noktasını ve bir uçuş süresi oranını kafadan kabaca hesapladım. Bu yüzden
ölçütleri sonuçtan değil, yalnız yayımlanmış ölçüm belirsizliklerinden kurdum. Uydurma (fit) yok; hiçbir değer sonuca göre
ayarlanmadı.

## Kısa sonuç

- **Denklemin itki halkası ölçümle doğrulandı.** F = γ·J_b·√(2mV/q), NASA NEXT itkicinin tezgâhta ölçülen 40 işletme noktasını
  serbest parametresiz öngörüyor:
  - 40 noktanın 40'ı 2σ bandında;
  - ortalama sapma +%1,0, kök ortalama kare %1,45;
  - 1000 V üstündeki 33 noktada kök ortalama kare %1,18.
- **Atom düzeyindeki halkalar ölçülmüş deneylerde doğrulandı:**
  - istek üzerine tek iyon (88,5 %, algılama sınırlı);
  - uçuş süresinin √m ölçeklemesi (1,97σ, sınırda);
  - yapay atomda I = e·f (0,27 ppm).
- **Programlanabilir atom → vakum yayımı halkası ölçülmüş bir aygıtta var.** SiOx'ta gerilimle taşınan oksijen, ns içinde açılıp
  kapanan bir elektron yayıcı kuruyor; 100'lük dizi 24 saat kararlı.
- **Ölçülmemiş tek halka hücreler** (SKY130, 4 K). İtki denklemi doğrulandığı için bu halka itki tezgâhı gerektirmiyor; 4 K'de
  akım basamaklarını ölçmek yetiyor.

## Sonuç tablosu

| no | halka | ölçülmüş veri | sonuç |
|---|---|---|---|
| MC1 | F = γ J_b √(2mV/q) + F_nötr | NEXT, 40 nokta (Diamant 2010 tezgâh itkisi; Pollard 2010 plüm γ) | ✓: 40/40 ≤ 2σ; ortalama +%1,01 (≤ %1,5); kök ortalama kare %1,45 (≤ %2) |
| MC2 | v = √(2qV/m) → t ∝ √m | Schnitzler 2009: Ca⁺ ~12,8 µs, CaO⁺ ~15,3 µs | ✓ sınırda: öngörü 1,1833, ölçüm 1,1953 ± 0,0061 (+%1,0; 1,97σ) |
| MC3 | atom atom, istek üzerine | Schnitzler 2009: istenen 1/2/4 iyon → 1/2/4 olay; 123/139 başarı | ✓: %88,5 başarı, algılama verimi 87(+7/−11)% ile uyumlu |
| MC4 | I = e·f (yapay atom) | Zhao 2017: Si MOS kuantum noktası, 1 GHz, ~160 pA | ✓: e·f = 160,2 pA; 0,27 ppm içinde |
| MC5 | programlanabilir atom vakum yayımını anahtarlar | Wu 2018: elektroform edilmiş SiOx nano-diyot | ✓: yayım yalnız atom yapısı değişince; açma 50 ns, kapama 130 ns |
| MC6 | hücreler (N·I_c) | SKY130 ölçümü yok | ölçülmedi: model düzeyi (önceden kabul) |

## MC1: NEXT itkisi, 40 nokta

**Yöntem.** Öngörü F = γ·J_b·√(2 m_Xe V_bps/e) + 0,41 mN.
- J_b, V_bps ve ölçülen itki: Diamant 2010 Tablo 2. Ölçüm belirsizliği %0,7–1,2.
- γ = 0,945 (TL01–22) ve 0,955 (TL23–40): Pollard 2010'un ExB ve Faraday sondasıyla ölçtüğü aralıkların (0,94–0,95 ve 0,95–0,96)
  ortaları.
- 0,41 mN: Diamant'ın ayrıca ölçtüğü nötrleştirici itkisi.

**Belirsizlik bütçesi** (her nokta için σ):
- tezgâh belirsizliği;
- γ için %1;
- V_b ≠ V_bps için 20 V (gerçek ışın gerilimi deşarj ve kuplaj gerilimleri kadar farklı; TT10 değerleri elimde yok);
- ana deşarjın nötr itkisi için 0,3 mN (yalnız şekilde verilmiş).

**Sonuç** (`closure.json`, ölçütler önceden kayıtlı):

| ölçüt | değer | eşik |
|---|---|---|
| (a) 2σ içindeki noktalar | 40/40 (%100) | ≥ %90 |
| (b) ortalama sapma | +%1,01 | ≤ %1,5 |
| (c) kök ortalama kare sapma | %1,45 | ≤ %2 |
| tanı: V_bps ≥ 1000 V'taki 33 nokta, kök ortalama kare | %1,18 | — |
| tanı: en büyük sapma | +%4,5 (TL02, 300 V; 2σ = %7,6) | — |

**Dürüst artık.**
- **Sistematik kayma.** Öngörü ortalama ~%1 yüksek.
- **Gerilimle büyüyen sapma.** Aynı J_b'de sapma gerilimle artıyor: 1,2 A'da 679 V'ta +%0,1 → 1800 V'ta +%2,4; 3,52 A'da
  1179 V'ta −%0,8 → 1800 V'ta +%0,5.
  - Pollard'ın ölçtüğü eğilimle aynı yönde: ışın sapması çarpanı β, V_b arttıkça kötüleşiyor.
  - Ben γ'yı kullanım başına sabit aldım; bu yüzden eğilim artığa geçiyor.
- **Tezgâhtan geri çıkan γ ile plüm ölçümü arasındaki fark.**
  - Tezgâh itkisinden geri çıkan γ 0,904–0,963.
  - J_b ≤ 2 A'de TL06–22'nin çoğunda (TL18–19 hariç) plüm bandının (0,94–0,95) alt kenarından %0,3–1,9 aşağıda; TL23–40'ta
    yüksek gerilimde bandın (0,95–0,96) %0,5–0,9 altında.
  - Diamant da farkların β yapısını izlediğini ve β için daha fazla plüm ölçümü gerektiğini yazıyor.
  - Diamant'ın kendi karşılaştırması (TT10'un gerçek V_b'si ve plümle düzeltilmiş γ) en çok %1,2 farkla kapanıyor (37/40 noktada
    ≤ %1,0).
- **Atom başına tek yük (tutarlılık, bağımsız sınama değil).** J_b·m/(e·ṁ_deşarj) her akım grubunda 0,872 / 0,933 / 0,934 /
  0,932 / 0,893–0,894 çıkıyor; Pollard'ın kullanım değerleriyle (%87 / 93 / 89) uyumlu. Bu, akımla atom akışı arasındaki "atom
  başına bir yük" ilişkisini doğruluyor; TT10 akışları bu kullanıma göre ayarlandığı için bağımsız bir sınama sayılmaz.

## MC2–MC5: atom düzeyi

- **MC2, uçuş süresi.** Aynı tuzaktan aynı alanla çıkarılan Ca⁺ ve CaO⁺ için denklem t(CaO⁺)/t(Ca⁺) = √(m_CaO⁺/m_Ca⁺) = 1,1833
  öngörüyor. Ölçülen oran 15,3/12,8 = 1,1953 ± 0,0061 (okuma ±0,05 µs). Fark +%1,0, 1,97σ: ölçüt sınırda geçti. Ek sistematik:
  iyonların kristaldeki başlangıç konumları farklı.
- **MC3, istek üzerine atom.**
  - İstenen 1, 2, 4 iyon algılayıcıda 1, 2, 4 olay veriyor.
  - 139 atışta 123 başarı (%88,5); ölçülen algılama verimi 87(+7/−11)%, yani kayıp çıkarımda değil algılamada.
  - Ölçülen 19,47 km/s hızla iyon başına itki biti 1,292·10⁻²¹ N·s, yayılımı 3,2·10⁻⁴ (6,3 m/s). Kinetik enerji 78,5 eV; makalenin
    "80 eV"u buradan geliyor, bağımsız bir girdi değil.
- **MC4, yapay atom turnikesi.** Silisyum MOS kuantum noktası pompası 1 GHz'te ~160 pA taşıyor (e·f = 160,2 pA). Akım, birincil
  gerilim ve direnç standartlarına izlenebilir ölçümde 0,27 ppm belirsizlik içinde e·f'ye eşit.
- **MC5, programlanabilir atom yayıcısı.** Grafen ya da metal elektrotlu nano aralıkta SiOx gerilimle elektroform ediliyor: SiOx → Si
  indirgenmesi, yani oksijen alanla taşınıyor ve aralıkta bir tünel diyodu kuruluyor.
  - 6,9 V'taki anahtarlamayla yayım başlıyor: 21,3 V'ta 3,2 µA (iletim 103,8 µA, %3,1); 40 V'ta verim %16,6.
  - Açma ~50 ns, kapama ~130 ns; anahtarlama sürelerine uyuyor.
  - Denetimler: elektroform edilmemiş aralıkta yayım yok; SiNx üzerinde 210 V'a kadar yok.
  - 100 diyotluk dizi: 24,8 V'ta 73,4 µA / (82 × 18 µm²) = 4,97 A/cm²; 5·10⁻⁶ Pa'da 24 saatte %3'ten az bozulma.
  - Üç eklemli kaskadla verim %47,6 (Li 2020).
  - **Sınır:** programlama ikili (koptu ya da oluşmadı); yayım başlangıçtaki direnç durumuna bağlı değil. Yayılan elektron, iyon
    değil. Exp24'ün O/HfO₂ sayıları (α = 0,507 eV/(V/nm)) bu ölçümle sınanmadı; doğrulanan, sınıfın kendisi (oksitte gerilimle
    taşınan oksijen vakum yayımını açıyor).

## Yol notu: her hücreye çip üstü programlanabilir nötrleştirici

İyon ışını aynı akımda elektronla nötrleştirilmeli. MC5'in ölçülmüş sayılarıyla (dizi 4,97 A/cm²; verim 21,3 V'ta %3,1, 40 V'ta
%16,6):

| ölçek | yayıcı alanı | güç | ışın gücüne oran |
|---|---|---|---|
| hücre (Exp23, 0,268 mA, 5 kV In⁺, ışın 1,34 W) | 5,4·10⁻⁵ cm² (73 × 73 µm) | 0,065–0,18 W | %4,8–13,7 |
| NEXT TL40 (3,52 A, 1800 V) | 0,71 cm² | 0,85–2,4 kW | %13–38 |

Hücre ölçeğinde her hücreye kendi çip üstü nötrleştiricisi sığıyor ve gücü makul. NEXT ölçeğinde birkaç on watt harcayan hollow-katoda
göre kazanç değil; yalnız 2,5–4 sccm (0,24–0,39 mg/s) nötrleştirici yakıtını kurtarır.

## Ne doğrulandı, ne kaldı

**Denklem (her terimin ölçülmüş dayanağıyla):**

  F = γ · N(t) · I_c · √(2 m V / q)    ve    I = e·f (metrelenmiş tek yük)

| terim | ölçülmüş dayanak |
|---|---|
| √(2mV/q): pozitif enerji qV'nin atom momentumuna dönüşmesi | MC1, MC2 |
| γ (ışın sapması ve çift yük) | Pollard 2010 plüm ölçümü; MC1'de tezgâhla %1–2 içinde |
| akım ↔ atom akışı, atom başına bir yük | MC1 tutarlılık, MC3 |
| tek yük metrelemesi I = e·f | MC4 |
| yayımı açıp kapayan programlanabilir atom | MC5 |
| N(t)·I_c, hücrelerin kuantumlu akımı | **ölçülmedi** (Exp17–23 model) |

**Exp28 düzeltmesi.** Exp28'de AR1 için "NSTAR η = 0,61 ✓" yazmıştım. Bu η, F/P = 2η/v_e tanımından hesaplandığı için bir özdeşlikti,
sınama değildi. Gerçek ölçüm sınaması bu deneydeki MC1.

**Kalan tek ölçüm** (laboratuvar gerekiyor; burada yapılamaz). İtki J ve V'den MC1 doğruluğuyla (~%1,4) çıktığı için hücre halkası
yalnız akım ölçümüyle sınanabilir; itki tezgâhı ve vakum odası gerekmez.
- **Düzen:** SKY130 hücre test yapıları, 4 K kriyojenik prob istasyonunda.
- **Kayıtlı öngörüler:**
  - hücre başına akım basamağı I_c = 0,166 mA (3,3 V) / 0,268 mA (5 V);
  - V_DD adımlarında k = 1, 2, 4 hücre (Exp19); zorla düşüşte 3 → 7 (Exp20);
  - +30 mV tohumla doğum yeri 10/10 (Exp18);
  - 300 K'de bölünme yok (Exp21).
- **Ölçülmüş benzer:** n-GaAs'ta 4,2 K'de safsızlık kırılımı filamanları görüntülenmiş (Brandl 1989). Silisyum SKY130 için ölçüm
  yok.

## Kaynaklar (bu oturumda tam metin okundu)

- Diamant, Pollard, Crofton, Patterson, Soulas, *Thrust Stand Characterization of the NASA Evolutionary Xenon Thruster (NEXT)*,
  NASA/TM-2010-216895 (AIAA-2010-6701), Tablo 2 ve denklem 18:
  https://ntrs.nasa.gov/api/citations/20100042402/downloads/20100042402.pdf
- Pollard, Diamant, Crofton, Patterson, Soulas, *Spatially-Resolved Beam Current and Charge-State Distributions for the NEXT Ion
  Engine*, NASA/TM-2010-216790, §3.7:
  https://ntrs.nasa.gov/api/citations/20100033115/downloads/20100033115.pdf
- Schnitzler ve ark., *Deterministic Ultracold Ion Source targeting the Heisenberg Limit*, PRL 102, 070501 (2009);
  arXiv 0902.0456.
- Zhao ve ark., *Thermal-error regime in high-accuracy gigahertz single-electron pumping* (2017);
  arXiv 1703.04795.
- Wu ve ark., *Silicon Oxide Electron-Emitting Nanodiodes* (2018); arXiv 1805.01602.
- Li ve ark., *A Cascade Electron Source Based on Series Horizontal Tunneling Junctions* (2020); arXiv 2009.02552.
- Brandl, Völcker, Prettl, Appl. Phys. Lett. 55, 238 (1989): 4,2 K'de n-GaAs akım filamanları.
