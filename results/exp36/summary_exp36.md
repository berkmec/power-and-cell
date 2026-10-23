# Exp36 — Fiziksel yüzey çekirdeği, arka-plan cezalı aktüatör ve tek-atış okumayla tam atom kontrolü

Tarih: 2026-10-17. Kod: `experiments/exp36_physical_kernel/`
- `uc36.py`: v1, 4×4 yüzey çekirdeği + skaler arka-plan cezası + iki onaylı tek-atış geri besleme;
- `uc36_v2.py`: v2, daha zayıf ceza + beş onaylı durma kuralı;
- `uc36_v3.py`: v3, tek değişkenli daha güçlü ceza;
- `uc36_v4.py`: v4, dondurulmuş yerel-alan güvenlik satırları;
- `uc36_v5.py`: v5, SEG6 (10 nm adım) ve SEG8 (7,5 nm adım) çözünürlük koşusu;
- `fig36.py`: `fig_exp36.png` istatistik figürü (4500×3000 px).

Öngörü kayıtları (her biri ilgili koşudan önce yazıldı; sonradan düzenlenmedi):
- `predictions_exp36.json` (RS1–RS6), `predictions_exp36_v2.json` (RV1–RV6), `predictions_exp36_v3.json` (RW1–RW6), `predictions_exp36_v4.json` (RX1–RX6), `predictions_exp36_v5.json` (RY1–RY6).

Sonuç dosyaları:
- `universal_control36.json` ve byte-özdeşi arşiv `universal_control36_v1.json`;
- `universal_control36_v2.json`, `universal_control36_v3.json`, `universal_control36_v4.json`, `universal_control36_v5.json`;
- `fig_exp36.png`: Pareto yolu, dilim adımına göre okuma kalitesi, desen sınıfına göre arka-plan bozulması ve SEG6 derinlik yanıtı.

## Soru ve tanımlar

Exp36, Exp35'in açık bıraktığı üç sınırı aynı görev setinde kapatmayı sınar:

1. **Fiziksel alan çekirdeği:** Hedef atom merkezli idealleştirilmiş Gauss kipleri yerine, z=0 yüzeyindeki dilim kapılarından gelen ve derinlikle zayıflayan alan kipleri kullanılır.
2. **Tek-atış okuma:** Kararlar analog ortalama doluluktan değil, sensör hatalı ikili anlık görüntülerden verilir. SS1 her turda %1 bit çevirme hatasıyla tek okuma, SS5 %5 hatayla üç okumanın çoğunluğunu kullanır; durmak için beş ardışık tahmini-tam desen gerekir.
3. **Arka-plan bütçesi:** Aktüatör çözümü hedef düzeltmelerle birlikte hedef dışı site enerjisi değişimini cezalandırır; başarı, toplam hedef dışı modal çevrilmenin kayıtlı %5 sınırını aşmamasını gerektirir.

Görevler Exp35 ile aynıdır: yedi 60 nm ada (tohum 0, 1, 2, 5, 6, 7, 8), adada kalibre adres sırasına göre en fazla 4 hedef atom ve MONO/ALT/RAND desenleri; toplam 21 görev, 72 görev-biti. Bit eşlemesi Exp24 ile aynıdır: bit 0 = B⁻ (n=1), bit 1 = B⁰ (n=0). Sıcaklık 4 K, geri besleme dilimi 600 tarama, son doğrulama taze RNG ile 1200 taramadır.

Aktüatör tabanı bir tekdüze gövde kipiyle karesel yüzey dilimleridir:

\[
B_{i0}=1,\qquad B_{iq}=\exp\left[-\frac{(x_i-x_q)^2+(y_i-y_q)^2}{2\sigma^2}\right]\exp(-z_i/50\,\mathrm{nm}).
\]

Her turda katsayılar

\[
\min_a \|B_Ta-d\|^2+\mu_{BG}\|B_{BG}a\|^2+10^{-3}\|a\|^2
\]

çözümünden gelir; dilim katsayıları ±0,80 V ile, etkin küresel V_B 0,805–1,006 V penceresiyle sınırlıdır. Bu, Poisson çözümlü bir yerleşim ya da SKY130 metal çizimi iddiası değildir; site-enerjisine eşdeğer volt bir elektrostatik vekil çekirdektir.

**Tam kontrol beyanı kuralı** bütün sürümlerde aynı kaldı: Bir çekirdeğin hem SS1 hem SS5 protokolü ALT+RAND'da marjinal ≥0,95, tam desen anlık oranı ≥0,65, modal başarı ≥%80 ve toplam hedef dışı çevrilme ≤%5 vermelidir. Analog okumalı ANA yalnız ablasyondur.

## Kısa sonuç

**Skaler ceza ya da sert güvenlik kısıtları tek başına yetmedi; yüzey dilimi çözünürlüğü belirleyici çıktı.**

- 4×4 SEG4 çekirdeğinde v2 hedef kalitesine ulaştı fakat arka planı %6,00/%6,14 ile aşamadı. Daha güçlü ceza (v3) ve aktif güvenlik (v4) arka planı indirdi ama hedef kalitesi çöktü.
- **SEG6 (6×6, 10 nm adım)** kayıtlı tam-kontrol kuralını geçti: SS1'de marjinal 0,982, tam desen 0,928, modal 13/14, toplam arka plan %4,70; SS5'te 0,981, 0,925, 13/14 ve %4,75.
- **SEG8 (8×8, 7,5 nm adım)** daha geniş marj verdi: SS1 ve SS5'in ikisi de ALT+RAND'da 14/14 modal, marjinal ≈0,999 ve toplam arka plan %3,95/%4,21.
- Sınır hâlâ görev bazında tekdüze değildir: SEG6'nın tek modal kaybı tohum 7 RAND'dır ve RAND sınıfı arka planı final çekirdeklerde %9,2–10,4 ortalamasında kalır; kayıtlı %5 ölçütü bütün görevlerin toplam ortalamasıdır, desen başına ayrı sınır değildir.

## Pareto yolu: neden v5 gerekti?

Aşağıdaki tabloda marjinal/tam desen/modal değerleri ALT+RAND (14 görev), arka plan bütün 21 görevin ortalamasıdır.

| Sürüm / protokol | Marjinal | Tam desen | Modal | Toplam arka plan |
|---|---:|---:|---:|---:|
| v1 SEG4-SS1 (μ=0,020; 2 onay) | 0,769 | 0,413 | 6/14 | %3,34 |
| v1 SEG4-SS5 | 0,795 | 0,392 | 6/14 | %3,30 |
| v2 SEG4-SS1 (μ=0,005; 5 onay) | 0,960 | 0,913 | 13/14 | %6,00 |
| v2 SEG4-SS5 | 0,943 | 0,843 | 12/14 | %6,14 |
| v3 SEG4-SS1 (μ=0,010) | 0,814 | 0,554 | 8/14 | %4,40 |
| v3 SEG4-SS5 | 0,885 | 0,698 | 10/14 | %4,48 |
| v4 SEG4-SS1 (aktif güvenlik) | 0,777 | 0,405 | 6/14 | %2,22 |
| v4 SEG4-SS5 | 0,768 | 0,377 | 6/14 | %3,31 |
| **v5 SEG6-SS1 (10 nm)** | **0,982** | **0,928** | **13/14** | **%4,70** |
| **v5 SEG6-SS5** | **0,981** | **0,925** | **13/14** | **%4,75** |
| **v5 SEG8-SS1 (7,5 nm)** | **0,999** | **0,996** | **14/14** | **%3,95** |
| **v5 SEG8-SS5** | **0,999** | **0,997** | **14/14** | **%4,21** |

Tanı özellikle tohum 7 ALT görevinde keskindir: SEG4 v2, iki bit-1 hedefinde istenen +220 meV yerine +144/+151 meV verebilir ve desen modal olarak yazılamaz; cezasız çözüm +204/+211 meV ile yazar ama arka planı %23,2 bozar. Deterministik ilk-çözüm taramasında SEG4'ten SEG6/SEG8'e hedef artığı 12,38 → 6,80/4,81 meV ve arka-plan RMS'i 44,01 → 34,01/27,96 meV olur. V5 bu doğrusal tanıyı Monte Carlo'da da doğrular.

## Final protokol sonuçları

Aşağıdaki kalite sütunları ALT+RAND ortalamasıdır; arka plan sütunu kayıtlı ölçütteki gibi bütün 21 görevin ortalamasıdır. SEG4 satırları arşivlenmiş v2 referansıdır, v5 içinde yeniden koşulmadı.

| Protokol | Marjinal | Tam desen | Modal (ALT+RAND) | Toplam arka plan |
|---|---:|---:|---:|---:|
| G (tek küresel V_B) | 0,630 | 0,000 | 0/14 | tanım gereği 0 |
| SEG4-ANA (v2) | 0,962 | 0,920 | 13/14 | %5,36 |
| SEG4-SS1 (v2) | 0,960 | 0,913 | 13/14 | %6,00 |
| SEG4-SS5 (v2) | 0,943 | 0,843 | 12/14 | %6,14 |
| SEG4-SS1-U (v2) | 0,920 | 0,848 | 12/14 | %15,43 |
| SEG6-ANA | 1,000 | 1,000 | 14/14 | %4,58 |
| **SEG6-SS1** | **0,982** | **0,928** | **13/14** | **%4,70** |
| **SEG6-SS5** | **0,981** | **0,925** | **13/14** | **%4,75** |
| SEG6-SS1-U | 0,944 | 0,921 | 13/14 | %8,21 |
| SEG8-ANA | 1,000 | 1,000 | 14/14 | %4,21 |
| **SEG8-SS1** | **0,999** | **0,996** | **14/14** | **%3,95** |
| **SEG8-SS5** | **0,999** | **0,997** | **14/14** | **%4,21** |
| SEG8-SS1-U | 0,941 | 0,909 | 13/14 | %7,12 |

SEG6 ve SEG8 ANA koşularında bütün 21 görev modal doğrudur. SEG8'in iki gürültülü protokolü de 21/21 modal başarıya ulaşır. SEG6'da iki gürültülü protokolün ortak tek kaybı tohum 7 RAND'dır (SS1 19, SS5 18 güncellemede durma onayı alamaz; yanlış durma değildir). Cezasız U koşularının arka planı %8,21/%7,12 olduğundan daha ince çekirdek bile arka-plan cezasını gereksiz kılmaz.

Geri besleme çabası düşüktür: SEG6-SS1/SS5 medyan güncellemesi 0,5/0,0, SEG8'de 0,0/0,0; bütün final protokollerinde yanlış durma 0'dır. SEG6 dilim katsayılarının %98,8'i, SEG8'de %99,7'si ±0,50 V içindedir; ±0,80 V klipsinde katsayı yoktur (en büyük mutlak değerler 0,757/0,766 V). Bütün etkin küresel gerilimler hücre penceresindedir.

Desen kırılımı izolasyonun kalan zayıf yönünü gösterir. Final SS1/SS5 arka planı MONO'da yaklaşık %0,3–1,5, ALT'de %2,4–3,0, RAND'da %9,2–10,4 aralığındadır. Toplam ölçüt MONO+ALT ortalamasıyla %5'in altında kalır; RAND'a ayrı %5 sınırı konsaydı hiçbir final protokolü geçemezdi. En yüksek tek görev bozulmaları tohum 1 RAND (%23,1 SEG6-SS1; %18,7 SEG8-SS1) ve tohum 7 RAND (%15,9/%18,8) görevlerindedir.

## Hükümler

### v1–v4: SEG4 Pareto sınırı

- **RS1–RS5 ✗, RS6 ✓ (v1).** μ=0,020 ve iki onaylı durma arka planı (%3,34 SS1) korur fakat hedef kalitesi 0,769/0,413 ve modal 6/14'te kalır.
- **RV1–RV3 ✗, RV4–RV6 ✓ (v2).** μ=0,005 ve beş onay hedefi çözer; üç ana protokol de yalnız toplam arka-plan sınırında kalır (ANA %5,36; SS1 %6,00; SS5 %6,14).
- **RW1–RW5 ✗, RW6 ✓ (v3).** μ=0,010 arka planı %4,40/%4,48'e indirir ama SS1 marjinali 0,814 ve modal 8/14'e düşer.
- **RX1–RX5 ✗, RX6 ✓ (v4).** Sert aktif güvenlik SS1 arka planını %2,22'ye indirir fakat hedef marjinali 0,777 ve modal 6/14 olur. Sorun tek skaler ödünleşim parametresiyle kapanmaz.

### v5 (RY1–RY6)

- **RY1 ✓.** SEG6-ANA marjinal 0,99984, tam desen 0,99955, modal 14/14 ve arka plan %4,58 verir.
- **RY2 ✓.** SEG6-SS1 marjinal 0,982, tam desen 0,928, modal 13/14, arka plan %4,70; medyan güncelleme 0,5 ve yanlış durma 0'dır.
- **RY3 ✓.** SEG6-SS5 marjinal 0,981, tam desen 0,925, modal 13/14 ve arka plan %4,75 verir; SS1 ile modal fark 0 puandır. Böylece adım beyanı için ek marjinal ≥0,95 koşulu da tutar.
- **RY4 ✓.** SEG8-SS1/SS5 ikisi de tam-kontrol bütçelerini geçer; SEG8-SS1 modal başarısı SEG6'dan yüksek (14/14'e karşı 13/14) ve arka planı daha düşüktür (%3,95'e karşı %4,70).
- **RY5 ✗ (kayıtlı biçimiyle).** SEG6'da sığ/derin ilk-tur marjinalleri %98,9/%96,0 (fark 3,0 puan), son SS1 değerleri %100/%98,1'dir (fark 1,9 puan). Son fark koşulu tutar; ilk farkın ≥10 puan olması koşulu, ince çekirdeğin derin hedefleri baştan daha iyi sürmesi nedeniyle tutmaz.
- **RY6 ✓.** İki çekirdekte de katsayıların en az %98'i ±0,50 V içinde, klips oranı %0 ve bütün etkin V_B değerleri pencerededir.

**Nihai hüküm:** Kayıtlı kuralla en kaba geçen çekirdek **SEG6 (10 nm adım)**; SEG8 daha geniş güvenlik marjı verir. Bu, Exp35'in analog okumalı ve hedef-merkezli idealleştirmesine göre fiziksel çekirdek + gürültülü tek-atış okuma + toplam %5 arka-plan bütçesinin aynı anda kapatılmasıdır.

## Yöntem sınırları

1. **Çekirdek hâlâ vekildir.** Gauss yanal profili ve exp(-z/50 nm) zayıflaması; gerçek dilim geometrisi, oksit/gövde Poisson çözümü, kapasitif bölünme, süreç değişkenliği ve SKY130 çizim kurallarıyla doğrulanmadı. 10/7,5 nm adım, bu modelde gerekli çözünürlük ölçeğidir; SKY130'da doğrudan üretilebilirlik iddiası değildir.
2. **Okuma sensörü idealleştirilmiştir.** Fiziksel anlık görüntüler Monte Carlo durumundan alınır ve bağımsız bit çevirme uygulanır. Ölçüm geri-etkisi, akım basamağı eşiği, bant genişliği/gecikme, kalibrasyon kayması ve okuma sırasında durum bozulması modellenmez. Beş onay kuralı bu sensör modelinde yanlış durmayı 0'da tutar; laboratuvar dedektöründe ayrıca sınanmalıdır.
3. **Arka-plan ölçütü toplam ortalamadır.** Her görevde hedef dışı modal çevrilme, o desenin G çalışma noktasına göre ve görev ortalamasıyla hesaplanır. RAND sınıfı ve bazı tek görevler %5'in belirgin üstündedir; desen başına veya görev başına izolasyon daha sert bir sonraki ölçüttür.
4. **Kayıt küçüktür.** Yedi ada, en fazla 4 hedef ve önceden nitelikli adresler kullanılır. Derinlik istatistiğindeki 72 görev-biti aynı atomların üç desende tekrarıdır, bağımsız 72 atom değildir. Daha büyük kayıt, adres seçimi ve tutma süresi ayrıca ölçeklenmelidir.
5. **Fiziksel laboratuvar kontrolü yapılmadı.** Sonuç, Exp24 Coulomb-camı + Exp35 örnekleyici ailesi içinde model-içi tam kontroldür. SKY130 hücresinin gerçek yerel alan ve tek-atış okuma halkaları ölçülmedi.

## Sıradaki iş

1. **Desen başına izolasyon:** RAND için ayrı ≤%5 ve tek görevde üst sınır koyan denetleyici; gerekirse sıralı yazma/zamanlama ve hedef etkileşim grafiğine göre alt desenlere bölme.
2. **Gerçek elektrostatik çekirdek:** Dilim geometrisini Poisson/gövde modeliyle ve ölçülebilir çapraz konuşma katsayılarıyla değiştirme; 10 nm adımın hangi fiziksel kapı yığınına karşılık gelebileceğini kapatma.
3. **Dedektör modeli:** Akım basamağı amplitüdü, eşik, örnekleme süresi, geri-etki ve kayma ile tek-atış hata bütçesini fiziksel okuma zincirine bağlama.
4. **Ölçek ve tutma:** Daha çok hedef/ada, yazma sonrası sürüşü azaltıp saklama ve koşullu bağlama (NAND) görevlerini aynı hata bütçesinde sınama.
5. **Laboratuvar hattı:** Exp33 protokolündeki ileri-V_B ön ölçümü ve geniş/balastlı SKY130 koşusu; ardından hücre gövdesinde tek atom basamağı ve yerel alan çapraz konuşması ölçümü.

**Sonuç cümlesi:** Exp36, tam kontrolün yalnız daha iyi geri beslemeyle değil yeterli yüzey-alan çözünürlüğüyle geldiğini gösterir: 15 nm SEG4 Pareto sınırında kalırken 10 nm SEG6, gürültülü tek-atış okuma altında kayıtlı toplam hata bütçesini kapatır; 7,5 nm SEG8 bütün desenleri modal olarak yazar. Kalan ana sınır RAND görevlerindeki yerel arka-plan bozulması ve çekirdeğin gerçek aygıt elektrostatiğine bağlanmasıdır.

## Kaynaklar

Exp24 (`summary_exp24.md`, `small_bodies.json`, `ap24_atoms.py`), Exp35 (`summary_exp35.md`, `predictions_exp35*.json`, `universal_control35.json`, `universal_control35_v3.json`, `uc35.py`, `uc35_v3.py`) ve Exp36'nın beş öngörü/beş sonuç JSON dosyası ile `fig_exp36.png`.


