# Exp26 — Doğal silisyumda bütün komşulara erişim: alan kodlaması, ayrıştırma, enerji bütçesi ve kesin sınama

Tarih: 2026-10-16.

Kod `experiments/exp26_natural_register/` içinde. `nr26.py` komutları:
- `law`: yalnız aşırı ince anahtar (D0);
- `run`: D+ kipinde alan kodlaması ve gerçekleşebilirlik;
- `bus`: çekirdekler arası veri yolu ve Lie rankı;
- `pulse`: seçici darbe sınaması;
- `fixed`: koşudan sonra yapılan düzeltme;
- `measured`: yalnız ölçülmüş değerlerle doğruluk–erişim dengesi.

Ayrıca:
- `vmrev.py`: kesin çok-spinli MREV-16 sınaması, `survey` alt komutuyla sızıntı taraması;
- `nr26_fig.py`: şekil.

Öngörüler her koşudan önce yazıldı:

| dosya | kayıt saati |
|---|---|
| `predictions_nr26.json` | 08:47 |
| `predictions_nr26_measured.json` | 09:07 |
| `predictions_nr26_vmrev.json` | 09:55 |
| `predictions_nr26_survey.json` | 10:01 |

## Düzeltme (10:05; §3 ve §6'nın sonuçlarını geri alıyor)

§3 ("t_c ≤ 0,40 ms ile %99,995") ve §6 ("ölçülmüş değerlerle 6 nm'nin tamamı ~%92 doğrulukla") yalnız **adreslenebilirliği**
hesaplıyordu. Her ikisi de ölçülmüş toplu T2'nin (25 s) her çekirdeğe ve yazılmış her bite uygulandığını varsayıyordu.

**Kesin çok-spinli benzetim (§7) bu varsayımı çürüttü.** Gerçek ²⁹Si konumları ve dipol bağlarıyla, Ladd'ın sonlu darbeli
MREV-16'sı ölçülmüş döngüde (1,568 ms) çalıştırıldı. Sonuç:
- tek bir 17,5 s'lik yazma süresinde bile yazmaçtaki bitlerin yalnız ~%31'i korunuyor;
- sızıntının ana kaynağı en yakın komşu ²⁹Si çiftleri (²⁹Si'lerin %17'si, |J| 52–574 Hz) ve onlara bağlı çekirdekler;
- benzetim ölçülmüş bozunmadan iyimser (17,5 s'de eşevrelilik 0,86, ölçülen 0,50); gerçek durum daha kötü olmalı.

Bu yüzden ölçülmüş değerlerle doğal Si'de 6 nm'nin tamamı **yazılamaz**. 16 kat hızlı döngüyle bile (uzatma) en yakın komşu
çiftlerinin üyeleri bitlerini kaybediyor. §3 ve §6'daki tablolar adreslenebilirlik sonuçları olarak duruyor, sonuçları geçerli
değil.

## Kısa sonuç

- **Adresleme mümkün.** Yönü seçilebilen gradyanla her ²⁹Si bir yönde yalıtılabiliyor (§2).
- **Depolama ölçülmüş koşullarda çöküyor.** Yazma, yüzlerce uzun MREV penceresi gerektiriyor; bu sürede yazılmış bitler
  korunamıyor (§7).
- **Ölçülmüş değerlerle doğal Si'de kalan erişim.** D0 kipinde, ayrıştırmasız, aşırı ince donmuş çekirdekle bitler korunuyor.
  Ama tek tek yazılabilen çekirdek ~69 (6 nm'dekilerin %3'ü). "Yarıçaptaki herkes" ölçütüyle yalnız ~0,8 nm. Evrensel kapı
  (elektron üzerinden) en güçlü ~5 komşuyla (Exp25).
- **Yanıt:** bütün komşulara erişen evrensel programlanabilir atom, ölçülmüş değerlerle saflaştırılmış ²⁸Si'de mümkün
  (Exp25: 50 ppm'de 6 nm'deki ²⁹Si'lerin hepsi, yapılandırmaların %98'inde). Doğal Si'de değil; engel en yakın komşu
  ²⁹Si–²⁹Si çiftleri. Saflaştırma bu çiftleri pratikte yok ediyor: çekirdek başına olasılık %0,1'de %0,4, 50 ppm'de %0,02.
- **Doğrudan ölçüm (§8).** Silisyumda tek bir programlanabilir yapay atom, bağlı olduğu bütün spinli komşuları evrensel olarak
  yönetiyor: Thorvaldson 2025 (3P + elektron, bütün çiftlerde CZ ≥ %99,3) ve Edlbauer 2025 (4P/5P, bütün çekirdekler, ²⁹Si sinyali
  yok). İkisi de saflaştırılmış ²⁸Si'de, yani Exp25–26'nın gösterdiği malzeme koşulunda.
- **Ölçümle doğrulama (§8).** Yayımlanmış deneyler ayrıca şunları gösteriyor:
  - tek bir atomun elektrik alanla programlanması (Asaad 2020);
  - donör elektronuyla tek ²⁹Si kontrolü (Pla 2014);
  - elektron + 2 çekirdek evrensel kapı seti (%99,95 / %99,37; Mądzik 2022);
  - elmasta merkez spin + 9 çekirdek (Bradley 2019).

  50 ppm'de bir donörün çevresindeki bütün spinli atomlar ~2 tane, yani ölçülmüş yazmaçlar ölçeğinde. Aynı model yöntemi
  elmastaki ölçülmüş küme boyutlarını (27 ve 8) büyüklük düzeyinde veriyor.

## Soru

Exp25'in sonucu şuydu: doğal Si'de 6 nm içindeki ²⁹Si'lerin hepsi, merkez atom üzerinden ölçülmüş eşevrelilik içinde
kontrol edilemiyor (UN3 ✓). Bu bir doğa yasağı mı, yoksa elektron eşevreliliğine dayanan protokolün sınırı mı? Doğal Si
önemli, çünkü SKY130 yongaları doğal silisyumdan.

## Protokol

Her parça ölçülmüş bir sayıya bağlı; ölçülmüş aralığın dışına çıkan her yer "uzatma" diye işaretli.

- **P1, yazma kipi D+.** Hücre donörü iyonlaştırıyor (Exp24/25 yük kontrolü varsayımı). Elektron olmadığı için elektron T1/T2
  sınırı da aşırı ince kalabalığı da yok.
- **P2, MREV-16 homonükleer ayrıştırma.** Ladd 2005 ölçümü: doğal Si'de T2 = 25 s, t_c = 1,568 ms, T2 ∝ t_c⁻², kaymalar
  1/3 ölçekleniyor. Model: T2(t_c, ω) = min(tavan, 25 s · (1,568 ms / t_c)² · (120 Hz / max(ω, 120 Hz))²).
- **P3, yönü seçilebilen alan gradyanı.** f_k = γG (R_k · g). Seçici darbe τ = T2/10, W = 4/(sτ). Çekirdek k, bir g yönü onu
  ötekilerden en az d* = W/(γG) ayırıyorsa yazılabilir.
- **P4, veri yolu.** Dipol ağı, J = 4,748 Hz nm³ (1 − 3cos²θ)/r³.
- **Yasalar.**
  - L1 (dipol bölgesi): N = 0,63 √(K/W), K = 1,19·10⁵ Hz.
  - L2 (tek yön): P(yalıtık) ≈ exp(−d* S_k), S_k = Σ_j 1/r_kj.

### Kullanılan ölçülmüş değerler (bu oturumda kaynaktan okundu)

| değer | kaynak |
|---|---|
| Doğal Si tek kristal, MREV-16: ²⁹Si T2 = 25,0 ± 0,2 s, t_c = 1,568 ms, oda sıcaklığı, RF ~120 Hz kayık, π/2 = 15 µs | Ladd vd., PRB 71, 014401 (2005) |
| T2 ∝ t_c^−2,00±0,2; etkin Hamiltonyen −Σ(ω_j/3)I_z; polikristalde 1/f tavanı ~8 s; ²⁹Si T1 ≈ 5 saat | Ladd vd. 2005 |
| Ayrıştırmasız ²⁹Si Hahn T2 ≈ 5,6 ms (doğal Si) | Dementyev vd., cond-mat/0308584 |
| Tek donör elektronu T1 ≈ 6 s (1,5 T), T1⁻¹ ∝ B⁵ | Morello vd., Nature 467, 687 (2010) |
| Cihazlarda 3 T altında T1⁻¹ ∝ B; cihazdan cihaza 100 kat fark | Tenberg vd. 2019, arXiv 1812.06644 (özet) |
| P donörü Raman 1/T1 = C T⁹, C = 1·10⁻⁷ s⁻¹K⁻⁹ (2,6–6 K) → T1(4 K) ≈ 38 s | Castner 1963, Morse vd. 2015 Tablo I |
| T1 = 3·10³ s (1,2 K); Orbach ΔE = 126 K | Tyryshkin vd., cond-mat/0303006 |
| ²⁹Si temas Stark katsayıları, 6 kabuk; ters simetri eşleri 1 V/µm'de ~48 kHz ayrılıyor | Park vd., PRL 103, 106802 (2009) |
| Mikromıknatıs gradyanı 0,2–0,8 mT/nm (cihaz simülasyonu) | Kawakami 2014 ek; arXiv 2303.04960 |

## Öngörüler ve sonuç

| no | öngörü | sonuç |
|---|---|---|
| NR1 | L1 yasası sayımı ±%40 içinde tutar | ✗: sayım 4,5–7,8 kat fazla (KL temas kuyruğu, §1) |
| NR2 | Cihaz T1e'si ve MREV ile anahtar < %5 | ✗: %11,7 |
| NR3 | Yön seçilebilen gradyanla ≥ %99 yalıtım (d* ≥ 0,010 nm); tek sabit yön < %5 | ✓: 20/20; %3,5 |
| NR4 | Ölçülmüş MREV noktasında katı ölçütle kapsama < %1 | ✓: %0,68 |
| NR5 | Tam kapsama için t_c ≤ 0,4 / 0,1 / 0,03 ms (tavansız / 25 s / 8 s) | kısmen: 0,399 ✓; 0,112; 0,063 ✗ (yalnız adresleme) |
| NR6 | Veri yolu bağlı; Lie rankı tam | ✓ |
| NR7 | Sürüş ≥ 0,1 mW; atom başına ≥ 1 µJ | ✓: 0,23 mW; 1,5 mJ |
| NR8 | Ölçülmüş noktada katı ölçütle R_all 2–3 nm | ✓ (sınırda): 2,0 nm (yalnız adresleme) |
| NR9 | 6 nm'de ≥ %99 kapsama τ ≥ 0,55 T2 ve tersinme 0,90–0,95 ister | ✓: 0,7 T2, 0,922 (yalnız adresleme) |
| NR10 | κ = 2 ve tersinme ≥ 0,95'te R_all 3–5 nm | ✓ (sınırda): 5,0 nm (yalnız adresleme) |
| NV1 | En yakın komşusu olmayan hedef: yazma ≥ 0,99, en büyük çapraz ≤ 10⁻², depolama sızıntısı ≤ 10⁻² | ✗: 0,79; 0,047; 0,51 (kümede güçlü çift var) |
| NV2 | ≥ 250 Hz'lik çiftin üyesi 17,5 s'de > 0,1 sızar | ✓: 0,99 |
| NV3 | Dipol kapalıyken DANTE-MREV yazması tek spin modelini verir | kısmen: hedef 1,0000 (3/3); çapraz 2 kümede Bloch'un 2 katı içinde, üçüncüde 3,9 kat daha iyi |
| NV4 | 1,568 ms'de çift üyelerinin yarısından, çifte bağlıların dörtte birinden fazlası > 0,1 sızar; ötekiler %90'ında < 0,01 | kısmen: %85 ✓, %75 ✓; ötekilerin %55'i > 0,1 ✗ |
| NV5 | 0,098 ms'de çift üyesi olmayanların < %10'u > 0,01 sızar; bazı üyeler hâlâ > 0,1 | ✓: %0; üyelerin %35'i |
| NV6 | Ölçülmüş döngüde bir yazmada bitini koruyan pay < %85; §6 iddiası tutmaz | ✓: %31 |

## 1. Yalnız aşırı ince anahtar (D0)

20 doğal yapılandırma, 6 nm içinde ortalama 2109 ²⁹Si, B = 1 T ∥ [1,2,3], `law.json`:

| W | tek dal | iki dal |
|---|---|---|
| 0,1 Hz | 1927 (%91,4) | %90,2 |
| 1 Hz | 1390 (%65,9) | %64,5 |
| 10 Hz | 531 (%25,2) | %24,8 |
| 40 Hz (cihaz T1e sınırı) | 246 (%11,7) | %11,5 |

- **NR1 neden tutmadı.** L1'i dipol bölgesi için türettim. Ama 6 nm içinde KL temas kuyruğu baskın: 6 nm'de 10–25 kHz, dipol
  ~73 Hz. Sayımlar √W eğimini koruyor; etkin katsayı K_eff ≈ 5–7·10⁶ Hz.
- **Ayrıştırmasız (statik çoklu çizgiler ±91 Hz).** ~69 çekirdek (%3,3). `d0_frozen_reach.json` ek bilgileri:
  - donmuş çekirdek oranı %92,6: |J| > 1 Hz olan her çift |J|'den fazla ayrık, yani eş çevirme bastırılıyor;
  - yarıçaptaki herkesin yazılabildiği R_all medyanı 0,8 nm (~2 atom);
  - 2 nm içinde %36, 3 nm içinde %22 yazılabilir.
- **Ters simetri eşleri.** %2,3'ü modelde tam dejenere; hücre alanı ayırır (Park 2009).
- **Uyarı.** KL çürüme uzunlukları gerçekçi asimptotik değerlerden uzun; temas kuyruğu ~2,5 kat abartılı olabilir.

## 2. D+ kipinde alan kodlaması (adresleme)

`run.json`:
- 20.000 yön ve yerel aramayla her çekirdek bir yönde yalıtılabiliyor; en kötü b_k = 0,0154 nm, d*₉₉ = 0,0192 nm.
- Tek sabit genel yön 0,010 nm'de %3,5 yalıtıyor.
- L2 yasası doğrulandı:

  | d* | sayım | L2 |
  |---|---|---|
  | 0,001 nm | 0,661 | 0,655 |
  | 0,003 nm | 0,281 | 0,265 |

- Düzeltme: [1,2,3] örgü yönüydü (düzlemler çakışık); referans genel yönlerle yeniden hesaplandı.

## 3. Gerçekleşebilirlik: adresleme için döngü süresi (depolama yok sayıldı; bkz. Düzeltme)

| T2 tavanı | en uzun t_c | γG | T2 | kapsama | B1 | ort. güç |
|---|---|---|---|---|---|---|
| yok (t_c⁻²) | 0,399 ms | 20 Hz/nm | 387 s | %99,995 | 15,5 mT | 0,23 mW |
| 25 s | 0,112 ms | 282 Hz/nm | 24,5 s | %99,84 | 55 mT | 2,9 mW |
| 8 s | 0,063 ms | 891 Hz/nm | 7,8 s | %99,84 | 98 mT | 9,1 mW |
| ölçülmüş, katı ölçüt | 1,568 ms | 20 Hz/nm | 25 s | %0,68 | 3,9 mT | 15 µW |

## 4. Veri yolu ve evrensellik (Lie rankı)

`bus.json`:
- Dipol ağı |J| ≥ 0,2 Hz'de 20/20 bağlı; merkezdeki en güçlü bağlı çekirdekten en çok 4 adım.
- 4 çekirdekli kümelerde Lie rankı 10/10 tam, hem seçici hem global kontrolle.
- D0'da elektron çekirdeklerin %99,93'üne ≥ 1 Hz ile bağlı.
- Bu, ilke düzeyinde kontrol edilebilirlik demek. Uzun yazma sırasında depolamanın çökmesi (§7) pratik evrenselliği sınırlıyor.

## 5. Zaman ve enerji (adresleme tasarımı için)

- Tavansız noktada: 357 yön × 38,7 s ≈ 3,8 saat, atom başına 1,5 mJ.
- Ölçülmüş noktada: 514 yön × 17,5 s ≈ 2,5 saat, atom başına 63 µJ.
- Darbe sınaması: hedef tersinmesi 0,988 (T2 = 10τ); çapraz uyarma 4/τ'da 1,0·10⁻⁵, 2/τ'da 4,4·10⁻³.

## 6. Yalnız ölçülmüş değerlerle doğruluk–erişim dengesi (adresleme; depolama yok sayıldı)

`measured.json`, ölçülmüş nokta (T2 = 25 s, kayma ≤ 120 Hz, γG = 120 Hz / R):

| tersinme | κ = 2 (çapraz 4,4·10⁻³) ile R_all | κ = 4 (çapraz 10⁻⁵) ile R_all |
|---|---|---|
| 0,988 | 3,0 nm | 2,0 nm |
| 0,954 | 5,0 nm | 4,0 nm |
| 0,922 | 6,0 nm | 4,5 nm |
| 0,847 | 6,0 nm | 6,0 nm |

Bu tablo, yazılmış bitlerin yazma süresince korunduğunu varsayıyor. §7 bu varsayımı çürüttü.

## 7. Kesin çok-spinli sınama (vmrev.py)

**Yöntem.**
- Gerçek ²⁹Si konumlarından 7 spinlik kümeler (128 boyut), tam seküler dipol Hamiltonyeni.
- MREV-16: MREV-8 ve x fazları ters çevrilmiş MREV-8; Ladd'ın sonlu darbeleri (π/2 = 15 µs, t_c = 1,568 ms).
- Seçici yazma: döngü başına bir küçük global darbe (DANTE), Gauss alanı π, 17,5 s.

**Uygulamanın doğrulanması.**
- Kayma ve bağ yokken döngü birim matris (hata 7·10⁻¹⁶).
- Etkin alan −z yönünde, ölçek 0,343; Ladd'ın −1/3'üyle uyumlu.
- İki spinde 1 s sonunda değişim:

  | J | 10 Hz | 100 Hz | 365 Hz | 573 Hz |
  |---|---|---|---|---|
  | değişim olasılığı | 1·10⁻¹³ | 2·10⁻⁷ | 1,2·10⁻³ | 3·10⁻² |

**Üç hedef kümesi** (`vmrev.json`):

| küme | yazma (dipol açık) | depolama sızıntısı 17,5 s | dipol kapalı: yazma / en büyük çapraz |
|---|---|---|---|
| en yakın komşusu yok (kümede 574 Hz'lik çift var) | 0,79 | 0,51 | 1,0000 / 2,1·10⁻³ |
| 261 Hz'lik çiftin üyesi | 0,71 | 0,99 | 1,0000 / 2,3·10⁻³ |
| iki en yakın komşu (574 ve 365 Hz) | 0,35 | 0,71 | 1,0000 / 1,0·10⁻⁵ |

**Mekanizma: üç spin etkisi.**
- Hedef ile 574 Hz'lik çift birlikteyken 17,5 s'de sızıntı 0,47.
- Çift çıkarılınca 0,0008; hedef ile çiftin tek üyesi 0,0014.
- Kalan ayrıştırma hatası güçlü çifti karıştırıyor, çift de ona bağlı komşuyu sürüklüyor.

**Döngü süresi** (darbe genişliği döngüyle orantılı; `vmrev_tc_scan.json`; 17,5 s sızıntısı):

| küme | 1,568 ms | 0,784 | 0,392 | 0,196 | 0,098 ms |
|---|---|---|---|---|---|
| çifte bağlı hedef | 0,48 | 0,15 | 0,025 | 0,006 | 0,0015 |
| 261 Hz çift üyesi | 0,99 | 0,69 | 0,23 | 0,088 | 0,023 |
| iki en yakın komşu | 0,63 | 0,90 | 0,67 | 0,99 | 0,45 |

**Tarama** (60 hedef, 7 spinlik kümeler, 17,5 s depolama; `vmrev_survey.json`). Hücreler, 0,1'den fazla sızan oran:

| sınıf (yazmaçtaki payı) | 1,568 ms | 0,392 ms | 0,098 ms |
|---|---|---|---|
| en yakın komşu çiftinin üyesi (%17,3) | %85 | %80 | %35 |
| çifte ≥ 5 Hz bağlı (%41,8) | %75 | %0 | %0 |
| ötekiler (%41,0) | %55 | %0 | %0 |
| **bitini koruyan yazmaç payı** | **%31** | **%86** | **%94** |

**Gerçekçilik denetimi** (`vmrev_coherence.json`). Toplu eşevrelilik ölçümle karşılaştırıldı: benzetim 17,5 s'de 0,86 ve 25 s'de
0,83 veriyor; Ladd'ın ölçümü sırasıyla 0,50 ve 0,37. Benzetim iyimser; gerçek sızıntı en az bu kadar.

**Sonuç.**
- Bütün yazmacı yazmak yüzlerce 17,5 s'lik MREV penceresi istiyor. Ölçülmüş döngüde tek pencerede bile bitlerin ~%69'u
  bozuluyor.
- 16 kat hızlı döngüde (0,098 ms, B1 ~63 mT, uzatma) çift üyesi olmayanlar pencere başına ≤ 10⁻² sızıyor. Ama çift üyelerinin
  %35'i hâlâ > 0,1 sızıyor. 514 pencere birikince onlar da bozulur.
- Daha yüksek mertebeli diziler (BR-24 vb.) sınanmadı.
- Çiftleri ayrı ele alan bir yol (D0'da donmuş depolama ve çifte özgü kapılar) açık bir aday.

## 8. Ölçümle doğrulama: yayımlanmış deneyler ve modelin onlara uygulanması

Bu oturumda kaynaktan okunan ölçümler:

| ölçüm | sistem | neyi doğruluyor |
|---|---|---|
| Asaad vd., Nature 579, 205 (2020): tek ¹²³Sb çekirdeği (spin 7/2) yalnız elektrik alanla kontrol ediliyor; mekanizma dört-kutup etkileşiminin modülasyonu ve örgü gerinimi; evre kaybı 0,1 s | silisyum cihazı | tek atomun alanla programlanması (8 düzey) |
| Pla vd., PRL 113, 246801 (2014): doğal Si'de donör elektronuyla tek ²⁹Si; A = 2,205 MHz, 1,77 T; T2* = 2,4 ms, Hahn T2 = 6,3 ms; sınırlayan öteki ²⁹Si'lerle dipol etkileşimi | doğal Si | merkez atom komşu ²⁹Si'yi yönetiyor; doğal Si'de ²⁹Si banyosu belirleyici (§7 ile aynı mekanizma) |
| Mądzik vd., Nature 601, 348 (2022): elektron + iki ³¹P çekirdeği (A = 95 ve 9 MHz), 730 ppm ²⁹Si'li ²⁸Si, 1,33 T; tek kübit %99,95, CZ %99,37, GHZ %92,5 | saflaştırılmış Si | merkez elektron + çekirdekler için evrensel kapı seti |
| Bradley vd., PRX 9, 031045 (2019): NV elektronu + 9 çekirdek (8 ¹³C + ¹⁴N); 45 çiftin hepsi dolanık, 7 kübite kadar dolanıklık, 63 s eşevrelilik | elmas, ¹³C yoğunluğu 1,94 nm⁻³ (doğal Si'de ²⁹Si: 2,34 nm⁻³) | merkez spinin ~10 komşuyla evrensel etkileşimi |
| Abobeih vd., Nature 576, 411 (2019): tek NV ile 27 ¹³C'lik kümenin atom ölçeğinde görüntülenmesi | elmas | merkez spinin 27 komşuyla etkileşimi ölçülmüş |

Modelin bu ölçümlerle karşılaştırılması (`verify_literature.json`):

1. **Elmas.** Aynı yöntem (örgü + bolluk + dipol, |A| = C_d √(1+3cos²θ)/r³) şunu veriyor:
   - |A| ≥ 8 kHz: 27,5 ¹³C (p10–p90: 21–34); ölçülen 27 görüntülenen.
   - |A| ≥ 25 kHz: 8,8 ¹³C (5–13); ölçülen 8 kontrol edilen.

   Eşikleri ben varsaydım; bu yalnız büyüklük düzeyi uyumu.
2. **Doğal Si.** Ölçülen tek ²⁹Si'nin T2'si (6,3 ms) hacimli malzemedekiyle aynı ve dipol banyosundan geliyor. Bu, §7'nin
   "doğal Si'de ²⁹Si–²⁹Si etkileşimleri belirleyici" bulgusuyla tutarlı.
3. **KL modeli Pla'nın 2,205 MHz sitesini vermiyor.** En yakın kabuklar 1,745 ve 2,981 MHz; yani iç kabuk ayrıntısında model
   %20–35 hatalı. D0 sayımları temkinle okunmalı; D+ ve dimer sonuçları KL'ye dayanmıyor.
4. **Saflaştırılmış Si.** 50 ppm'de 6 nm içindeki ²⁹Si sayısı ortalama 2,2; donörlerin %61'inde ≤ 2, %92,5'inde ≤ 4.
   - Yani "etrafındaki bütün spinli atomlar" kümesi, evrensel kontrolü ölçümle gösterilmiş yazmaçlar ölçeğinde: Si'de
     elektron + 2 çekirdek, elmasta elektron + 9 çekirdek.
   - 730 ppm'de (Mądzik'in malzemesi) bu sayı 33, doğal Si'de 2118; bu ölçeklerde ölçüm yok.

**Ölçümle doğrulanan kısım:**
- Atomun alanla programlanması ölçülmüş (Asaad 2020).
- Merkez atomun çevresindeki spinli atomlarla evrensel etkileşimi ölçülmüş; ölçek Si'de ~3, elmasta ~10 spin.
- 50 ppm ²⁸Si'de bir donörün çevresindeki bütün spinli atomlar (~2) bu ölçekte. Bileşenler ayrı ayrı ölçülmüş: donör
  elektronuyla ²⁹Si kontrolü (Pla 2014) ve elektron + 2 çekirdek evrensel kapı seti (Mądzik 2022).
- Tam olarak "50 ppm'de bir donör + onun bütün ²⁹Si'leri" deneyi yapılmadı. Bu benim modelimin bileşenlerin ölçümlerinden
  çıkardığı sonuç.
- Doğal Si'de bütün komşulara erişim hiçbir ölçümde yok. Exp26 modeline göre ölçülmüş ayrıştırmayla da olamaz (dimerler).

**Doğrudan ölçüm: tek bir yapay atom bağlı olduğu bütün spinli komşuları evrensel olarak yönetiyor.** İki deney, kaynak metinden
okundu:

| ölçüm | yapay atom | "bütün komşular" | ölçülen kontrol |
|---|---|---|---|
| Thorvaldson vd., Nature Nanotechnology 2025 (arXiv 2404.08741) | 3 ³¹P çekirdeği + 1 elektron, 1,5 nm²'ye hassas yerleştirilmiş, ~200 ppm ²⁹Si'li ²⁸Si, 1,45 T, 15 mK | Hamiltonyende yalnız 3 ³¹P çekirdeği (A = 6, 68, 103 MHz); ²⁹Si bağı bildirilmemiş | tek kübit %99,94–99,98 (4 kübitin hepsi); **bütün** çekirdek çiftlerinde CZ: %99,65, %99,49, %99,32; GHZ (3 çekirdek) %96,2; Grover %93,5 (ortalama %89,4) |
| Edlbauer vd., Nature 2025 (arXiv 2506.03567) | 4P ve 5P kayıtları; her kayıtta çekirdekler ortak bir elektrona aşırı ince bağlı; iki elektron değiş-tokuşla bağlı (11 kübit) | iyonlaşmış NMR'da ²⁹Si'den sinyal yok; "bütün çekirdek spinleri için" tek ve çok kübit kontrol | tek kübit ≥ %99,9 (en çok %99,99); çekirdekler arası CZ %99,90; elektron–elektron CROT %99,64; GHZ (3 kübit) %92; 8 çekirdek spinine kadar dolanıklık |

- **Bunlar neyi doğruluyor.** Silisyumda tek bir programlanabilir yapay atom (fosfor kümesi ve onun elektronu), etkileştiği bütün
  spinli komşu atomları yönetebiliyor; evrensel kapı seti, her çift arasında iki kübit kapısı ve çok parçalı dolanıklıkla.
  İki yapay atom değiş-tokuşla birbirine bağlanabiliyor.
- **Exp25–26 ile tutarlılık.** Malzeme koşulu benim bulduğumla aynı: saflaştırılmış ²⁸Si (≤ 200 ppm ²⁹Si). Burada atomun spinli
  komşuları birkaç tane ve hepsi yönetilebiliyor. Doğal Si'de (%4,7) bu ölçekte bir sonuç yok; Exp26 nedenini ²⁹Si çiftleriyle
  açıklıyor.
- **Fark.** Buradaki "atom" tek P değil, P kümesi (yapay atom). Komşuları kümenin kendi çekirdekleri; spinsiz ²⁸Si'ler yük ve
  orbital yoluyla etkileşiyor, spin yoluyla programlanmıyor.

Ek ölçüm, Mądzik vd., Science Advances 6, eaba3442 (2020), saflaştırılmış ²⁸Si'de tek ³¹P donörü: donör elektronu yakındaki
²⁹Si'lerin hareketini donduruyor; elektron çıkarılınca hareket geri geliyor ve T2* olağan değerine iniyor. Bu, Exp26'nın
mekanizmasıyla aynı: D0'da donmuş çekirdek bitleri korur (%92,6 donmuş), D+'da çekirdek hareketi yazma sırasında bitleri bozar
(§7).

## 9. Ölçüm önerisi: tek donörün bütün komşularını yönettiğinin laboratuvarda gösterilmesi

Henüz hiçbir yayımlanmış ölçüm, tek bir programlanabilir atomun bütün komşularını evrensel olarak yönettiğini göstermiyor. Bu
çalışmanın öngörülerini doğrudan sınayacak deney:

1. **Cihaz.**
   - 50 ppm ²⁹Si'li ²⁸Si epikatmanda iyon ekimiyle tek ³¹P donörü;
   - tek elektron transistörüyle tek atış spin okuma (Morello 2010 tipi);
   - B0 ≈ 1–1,5 T, ~100 mK;
   - ESR anteni ve donöre yakın elektrik kapısı.
2. **Komşuların tamamını bulma.** Elektronun Ramsey ve dinamik ayrıştırma (CPMG/XY-8) spektroskopisiyle aşırı ince bağlı bütün
   ²⁹Si'ler aranır (elmasta Abobeih 2019 yöntemi).
   - Öngörü: 6 nm içinde ortalama 2,2 çekirdek; donörlerin %61'inde ≤ 2.
   - Algılama eşiği ≤ 1 kHz olmalı.
3. **Her komşu için tam kapı seti.**
   - elektronla tek atış çekirdek okuma (Pla 2014);
   - seçici NMR;
   - elektron üzerinden geometrik faz ile CZ (Mądzik 2022);
   - kapı kümesi tomografisi.
4. **"Bütün komşular" ölçütü.** Elektron + ³¹P + bulunan bütün ²⁹Si'lerin GHZ durumu ve tam süreç tomografisi. Karşılaştırma
   ölçütü Mądzik 2022'nin 3 kübitlik GHZ'si (%92,5).
5. **Doğal Si öngörüsünün sınanması.** Doğal Si ve %0,1 ²⁹Si numunelerinde toplu ²⁹Si NMR yapılır; MREV-16 t_c = 1,568 ms ve
   0,4 ms'de uygulanır.
   - Öngörü: doğal Si'nin eşevrelilik eğrisinde ~%17'lik hızlı sönen bir bileşen (en yakın komşu çiftleri).
   - %0,1 numunede bu bileşen olmamalı.
6. **D0/D+ donma testi.** Mądzik 2020 yöntemiyle, elektron varken ve yokken yakın ²⁹Si bitlerinin ömrü ölçülür. Exp26'nın D0
   depolama ve D+ sızıntı öngörüleriyle karşılaştırılır.

## Sınırlar

1. Model + literatür; benim tarafımdan ölçülmedi. Doğrulama yayımlanmış ölçümlerle karşılaştırma düzeyinde (§8). Asıl doğrulama
   §9'daki laboratuvar deneyi.
2. Çok-spinli benzetim 7 spinlik kümelerle yapıldı; daha büyük çevre sızıntıyı artırabilir. Ölçülmüş eşevrelilikle
   karşılaştırma, benzetimin iyimser olduğunu gösteriyor.
3. Başlatma ve okuma modellenmedi.
4. ²⁸Si ve ³⁰Si spinsiz; bu yolla programlanamaz.
5. D0 sayımları KL temasına dayanıyor; temas kuyruğu abartılı olabilir.

## Sonuç

- **Doğal silisyum, ölçülmüş değerlerle.**
  - Merkez atom bütün ~2100 ²⁹Si'yi tek tek adresleyebilir; buna gradyan kodlaması gerekiyor (§2).
  - Ama ölçülmüş en iyi ayrıştırma (MREV-16, 1,568 ms) yazma süresince bitleri koruyamıyor. Bunu kesin çok-spinli benzetim
    gösteriyor ve benzetim ölçümden iyimser.
  - Ölçülmüş değerlerle kalan: D0 kipinde ~69 çekirdek ve Exp25'teki elektron kapılarıyla en güçlü ~5 komşu.
  - Engel en yakın komşu ²⁹Si–²⁹Si çiftleri: ²⁹Si'lerin %17'si, |J| 52–574 Hz. 16 kat hızlı döngüyle bile tam çözülmüyor.
- **Saflaştırılmış silisyum.** Exp25: 50 ppm'de 6 nm içindeki bütün ²⁹Si'ler ölçülmüş T2e içinde, yapılandırmaların %98'inde.
  Çiftler pratikte yok; evrensel programlanabilir atom ölçülmüş değerlerle (model düzeyinde) mümkün.
- **Yanıt:** "bütün komşulara erişen evrensel atom" için doğru malzeme saflaştırılmış ²⁸Si. Doğal Si'de bu, ölçülmüş
  eşevrelilikle yapılamıyor; sebebi ²⁹Si dimerleri.
- **Ölçümle durum (§8).** Silisyumda tek bir programlanabilir yapay atomun bağlı olduğu bütün spinli komşuları evrensel olarak
  yönettiği ölçülmüş.
  - Thorvaldson 2025: 3P + elektron; bütün çiftlerde CZ ≥ %99,3; GHZ %96,2.
  - Edlbauer 2025: 4P ve 5P kayıtları; bütün çekirdekler kontrol ediliyor; CZ %99,90; ²⁹Si sinyali yok.

  Malzeme koşulu Exp25–26'nın bulduğuyla aynı: saflaştırılmış ²⁸Si. Bileşenleri de ölçülmüş:
  - alanla atom programlama (Asaad 2020);
  - donör elektronuyla ²⁹Si kontrolü (Pla 2014).

  Tek bir ³¹P'nin kendi ²⁹Si'lerinin hepsini yönettiği deney ise henüz yok (§9).
- **Sınanabilir öngörüler:**
  1. Doğal Si'de MREV-16 altında ²⁹Si eşevreliliği bileşenlere ayrılmalı: en yakın komşu çiftlerinden gelen ~%17'lik hızlı
     sönen bir bileşen bekleniyor.
  2. ≤ %0,1 ²⁹Si'de bu bileşen kaybolmalı.

## Dosyalar

- Öngörüler: `predictions_nr26*.json`.
- Sonuçlar:
  - `law.json`, `run.json`, `bus.json`, `pulse.json`, `measured.json`;
  - `vmrev.json`, `vmrev_tc_scan.json`, `vmrev_survey.json`, `vmrev_coherence.json`;
  - `d0_frozen_reach.json`, `verify_literature.json`.
- `fig_exp26.png`, `*.log`.
