# Exp25 — Atom, alan ve enerji manipülasyonu: evrensel programlanabilir atom

Tarih: 2026-10-16. Kod: `experiments/exp25_universal_atom/`
- `ap25_energy.py`: enerji;
- `ap25_logic.py`: alanla bağlı atom mantığı;
- `ap25_register.py`: evrensel merkez atom;
- `ap25_fig.py`: şekil.

Tahminler her koşudan önce yazıldı (`predictions_*.json`). Şekil: `fig_exp25.png`.

## Hedef ve yorum

Kullanıcı: bütün veriyi atom, alan ve enerji manipülasyonuna odakla; önce enerjiyle dene; sonra programlanabilir atomu
evrenselleştir, öyle ki etrafındaki bütün atomlarla etkileşsin.

Benim yorumum:
- **Enerji:**
  - hücre enerjisinin yönlendirilmesi (E1);
  - atom işlemlerinin enerji maliyeti (E2);
  - hücre alanıyla atomun enerji düzeylerinin kaydırılması (E3).
- **Evrensel atom, klasik:** her atom Coulomb alanıyla bütün komşularını etkiliyor; bundan evrensel kapı (NAND) kurulması.
- **Evrensel atom, kuantum:** bir merkez atomun (³¹P donör) aşırı ince bağla bütün komşu çekirdeklerle etkileşmesi ve onları
  tam kapı setiyle (ESR + seçici NMR + koşullu CNOT) yönetmesi.

## Kısa sonuç

1. **Enerjiyle doğrudan atom oynatmak işe yaramıyor, enerji düzeylerini oynatmak yarıyor.**
   - Hücre enerjiyi yoğunlaştırmıyor (V_M kenetlenmesi, Exp23).
   - Ölen hücrenin gücü %97–98 verimle kalanlara geçiyor, ama payı yakınlık belirliyor (%73/%27). Kapı ön gerilimi bu payı
     yalnızca ±%4 değiştiriyor.
   - Hücrenin alanı (0,4–1,8 V/µm) donör atomunun aşırı ince ve Zeeman enerjilerini kHz–MHz kaydırıyor. Bu, atomu frekansla
     seçmek için yeterli.
2. **Alanla bağlı evrensel mantık çalışıyor (model, 4 K).**
   - Bir atomun alanı 30 nm içindeki bütün atomlarda kT'yi aşıyor.
   - Üç doğal bor atomundan kurulan NAND kapısı, eşikler her giriş durumu için kalibre edilince 28 üçlünün 27'sinde, 10
     adanın 9'unda çalışıyor (dört satır ≥ 0,9).
   - İlk yöntemimde pencere hatası vardı: 3/28, LG1 ✗. Bu kayıtta duruyor.
3. **Kuantum evrensel atom çalışıyor (model).**
   - Merkezdeki ³¹P elektronu bütün yakın ²⁹Si çekirdekleriyle etkileşiyor.
   - 1 T ve 10 kHz seçicilikle doğal Si'de donör başına ~5 çekirdek tek tek yönetilebiliyor.
   - Tek RF bobiniyle elektron + 4 çekirdek GHZ dolanıklığı kuruluyor: doğruluk 0,999, öteki ~2000 çekirdek dahil edilince
     0,987. T2e ≥ 10 ms ile ~0,91.
   - Sınır: doğal Si'de güçlü bağlı komşuların yalnızca %1–3'ü tek tek seçilebiliyor (spektral kalabalık). Saflaştırma bu
     oranı %40–64'e çıkarıyor ama sayıyı düşürüyor.
   - Doğal Si'deki 100 kHz'lik ilk tasarım tutmadı (RG1, RG4 ✗).

## 1. Enerji manipülasyonu

**E1, enerji yönlendirme** (Exp18 3,3 V durumu, sağ hücre 30 s rampayla öldürüldü; `energy.json`):

| durum | uzak hücrenin payı | yakın hücrenin payı | verim | doğum/ölüm |
|---|---|---|---|---|
| yönlendirmesiz | 0,266 | 0,734 | 0,980 | 0 |
| uzağa +30 mV | 0,257 | 0,743 | 0,978 | 0 |
| yakına +30 mV | 0,222 | 0,778 | 0,973 | 0 |
| uzağa +15 mV | 0,231 | 0,769 | 0,979 | 0 |
| yakına +15 mV | 0,253 | 0,747 | 0,976 | 0 |

Tahminler:
- EN1 ✓ (iki hücre de ≥ %20 alıyor).
- EN2 ✗ (+30 mV ile ≥ %80 olmadı; uzağa yönlendirme hiç işlemedi).
- EN3 ✓ (verim 0,97–0,98; doğum ve ölüm yok).
- EN4 ✗ (uzak hücrede monoton değil).

Sonuç: enerjinin nereye gideceğini kapı değil geometri belirliyor. Yön seçmenin yolu, hücrelerin nerede doğacağını programlamak
(Exp18); bu yeni bir sınama ister, yapılmadı.

**E2, enerji maliyeti:**

| işlem | enerji |
|---|---|
| Landauer sınırı, 4 K | 0,24 meV |
| tek atomun yükünü çevirme (8 kT/e'lik basamak) | 2,8 meV (11,5 × Landauer) |
| iyonik yazma, 1 µs | ~2,7·10⁸ eV/bit (atom başına 1,4·10⁸ eV; bariyer ~1 eV) |
| hücreyi çalışır tutmak, 136 µW × 1 ms | ~8,5·10¹¹ eV |
| ³¹P çekirdek çevirme | foton 7·10⁻⁸ eV; 1 mW × 50 µs RF darbesi ~3·10¹¹ eV |

Atomun kendi ihtiyacı küçük; maliyeti sistem getiriyor.

**E3, enerji düzeyi manipülasyonu** (kuadratik Stark; Bradbury 2006, doğrulandı):
- Hücre cephesinde E = 0,9 V/genişlik.

| cephe genişliği | alan | ³¹P aşırı ince kayması (kuram, −2·10⁻² µm²/V²) | Sb benzeri (ölçülen, −3,7·10⁻³) | ESR (g-faktörü, 1 T) |
|---|---|---|---|---|
| 2,2 µm | 0,41 V/µm | 0,39 MHz | 0,07 MHz | 47 kHz |
| 1,0 µm | 0,90 V/µm | 1,9 MHz | 0,35 MHz | 227 kHz |
| 0,5 µm | 1,80 V/µm | 7,6 MHz | 1,4 MHz | 906 kHz |

- Alanlar iyonlaşma eşiğinin (3,7 V/µm) altında. Kaymalar kHz'lik çizgi genişliklerinden çok büyük; hücre, atomun rezonansını
  seçerek onu öteki atomlardan ayırabilir (Kane "A-kapısı").

## 2. Alanla bağlı atomlar: evrensel mantık (NAND)

`ap25_logic.py`; ilk tahminler `predictions_logic.json`, düzeltilmiş yöntem `predictions_logic2.json`.

**Kural:** çıkış atomu o, ancak V_b < V_th(n_i, n_j) iken iyonize (bit 1). Girişlerin yükü eşiği W_oi n_i + W_oj n_j kadar
kaydırıyor. V_b, V_th(1,1) ile (0,0), (0,1), (1,0) eşiklerinin en küçüğü arasında seçilirse çıkış = NAND(giriş₁, giriş₂).

**Model:**
- 60 nm'lik adalar (10 tane), 4 K Metropolis.
- Girişler yerel kapılarla ±80 mV zorlanıyor. Bu bir varsayım; SKY130 ~20 nm'lik kapı çizemez.
- Arka plandaki atomlar serbest.

**İlk deneme (LG1 ✗):** 28 üçlünün 3'ü, 10 adanın 3'ü. Başarısız 25 üçlünün 23'ünde yalnızca (1,1) → 0 satırı çöktü. Nedeni
benim yöntem hatam: pencereyi girişlerin referans eşikte nötr olduğu varsayımıyla kurmuşum.

**Düzeltilmiş yöntem** (çıkışın eşiği her giriş durumu için ayrı kalibre):
- LG4 ✓: 10 adanın 9'unda NAND var.
- LG5 ✓: penceresi ≥ 5,5 mV olan 25 üçlünün 25'i çalışıyor.
- Toplamda 28 üçlünün 27'si çalışıyor; pencereler 1,8–33 mV, satır doğrulukları çoğunlukla 0,96–1,00.

**Erim (LG3 ✓):** bir atomun alanı 30 nm içindeki atomların %100'ünde kT'yi aşıyor. Programlanabilir atom bütün komşularıyla
etkileşiyor.

**Anlamı:** NAND tek başına evrensel olduğundan, alanla bağlı doğal atom üçlüleri her klasik mantığı kurabilir. Sınırlar:
- yerel kapı varsayımı;
- kalibrasyon şart;
- model düzeyi; ölçüm yok.

## 3. Kuantum evrensel atom: merkez donör ve ²⁹Si komşuları

`ap25_register.py`; tahminler `predictions_register.json` ve `predictions_register2.json`.

**Model:**
- Kohn–Luttinger temel durumu: a = 2,5 nm, b = 1,4 nm, k₀ = 0,844·2π/a₀ (Koiller 2004).
- Temas aşırı ince bağı, ölçülen [100] ilk kabuk değerine (2,981 MHz, Hale ve Mieher 1969) kalibre. KL'nin doğruluğu ~%50 (Feher
  1959).
- 1 T'de çekirdek frekansı f_k(±) = |ν_n ± A_k/2|. Bir çekirdek, f_k(+) bütün öteki geçişlerden ≥ Δ_min uzaksa tek tek
  adreslenebilir.
- GHZ: tek RF bobini (bütün çekirdekler her darbeyi görüyor), tam dönen dalga yaklaşımı, yerel Z fazlarına kadar doğruluk.

**Doğal Si, 100 kHz (ilk tasarım):**
- RG1 ✗: ortalama 0,82 adreslenebilir çekirdek; öngörü ≥ 3.
- RG4 ✗: 2 nm içindeki çekirdeklerin hepsi yalnızca yapılanmaların %4'ünde > 1 kHz bağlı; vadi girişimi düğümleri var.
- RG2 ve RG3 sınanamadı (4 adreslenebilir çekirdekli yapılanma yok).

**Tasarım taraması** (`register_isotopes.json`; RG5–RG7 koşudan önce kaydedildi):

| ²⁹Si | Δ_min | adreslenebilir/donör | ≥ 4 olan | 3 nm içindeki güçlülerden adreslenebilen |
|---|---|---|---|---|
| %4,685 | 100 / 30 / 10 kHz | 0,66 / 2,34 / 5,12 | %0 / %18 / %84 | %1,2 / %1,8 / %2,7 |
| %1 | 100 / 30 / 10 kHz | 0,54 / 2,02 / 4,48 | %0 / %18 / %64 | %5 / %7 / %11 |
| %0,1 | 100 / 30 / 10 kHz | 0,30 / 1,06 / 2,36 | %0 / %2 / %14 | %40 / %40 / %48 |
| %0,05 | 100 / 30 / 10 kHz | 0,24 / 0,96 / 1,58 | %0 / %2 / %2 | %61 / %57 / %64 |

- RG5 ✗ (%0,1'de 1,06 çekirdek ve %40; öngörü ≥ 3 ve ≥ %80).
- RG6 ✓, küçük örneklem (n = 1 ve 7).
- GHZ doğruluğu (elektron + 4 çekirdek): eşevrelilik kaybı olmadan 0,99–0,999; T2e = 10 ms ile 0,905–0,966; T2e = 0,3 ms ile
  ~0,50–0,55.
- Seyirci düzeltmesi (kayıtlı ölçütün dışında, ayrıca hesaplandı): öteki bütün çekirdeklerin küçük dönmeleri dahil.
  - Saflaştırılmış Si'de 0,980–0,995.
  - Doğal Si'de, 10 kHz ve 12 yapılanmada 0,987 (en kötüsü 0,978).

**Hücrenin rolü:** hücre gerilim adresiyle merkez atomun yükünü seçiyor. D⁰'da atom bütün komşularına bağlı (hesaplama), D⁺'da
bağ sıfır (depolama). Stark kaydırması (E3) birden çok merkez atomdan hangisinin rezonansta olacağını seçiyor. Bunun için
n-gövdeli ayna hücre gerekiyor (varsayım).

## 4. Bütün komşular: merkez atom üzerinden dolaylı evrensel kontrol

`ap25_control.py`; tahminler `predictions_control.json` ve `predictions_control2.json`, koşudan önce yazıldı.

**Fikir (evrensel kuantum arayüzü, Lloyd–Landahl–Slotine 2004):** tek tek frekans seçimi spektral kalabalıkta tıkanıyor. Ama
yalnız global kontroller ve kayma Hamiltonyeni H₀ = Σ_k A_k S_z I_z^k ile üretilen Lie cebiri tam su(d) ise sistem tam kontrol
edilebilir. Global kontroller: elektron ESR x/y ve bütün çekirdekleri aynı anda süren tek RF x/y. Bunun için bütün bağların
farklı ve sıfırdan büyük olması yeter. Bedeli zaman: süre ~ 1/(en küçük bağ farkı).

**Bağlar:** temas kısmı (Kohn–Luttinger, kalibre) + seküler dipol kısmı T_k(3cos²θ_k − 1).

**Sonuçlar:**
- CT1, kısmen ✗:
  - 3 çekirdekli kümelerin 10/10'u tam (255/255).
  - 4 çekirdekli kümelerin 7/10'u tam; 3 kümede rank 638 ve 318/1023. Bu kümelerde bağlar birebir eşit; B ∥ [001] simetri
    eksenindeyken eşdeğer konumlar ayrışmıyor.
- CT2 ✓: birebir eşit çift → 38/63. Eşdeğer [100] konumları dipolle 3177 ve 2883 kHz'e ayrılıyor → 63/63.
- CT5 ✓: B ∥ [1,2,3] iken aynı 10 kümenin 10'u tam (1023/1023); ters simetri çifti yok.
- CT3 ✗ / CT6 ✓: 1 kHz aralıklı (seçilemeyen) çekirdeğe CNOT, T = 4 ms. 400 yinelemede 0,989 (eşik 0,99); 1500 yinelemede
  1,0000 ve 0,9989.
- CT4 ✓: elektron + 3 çekirdeğin hepsini kapsayan GHZ (ikisi 1 kHz aralıklı) 1,0000 doğrulukla hazırlanıyor.
- CT7 (keşif): doğal Si'de en güçlü 10 komşu arasındaki en küçük bağ farkı medyan 1,3 kHz (p10: 0, p90: 5,4 kHz), yani
  kontrol süresi ~0,8 ms. p10 = 0 değeri, bu gerçel KL modelindeki ters simetri eşitliklerinden geliyor; gerçek donörde ters
  simetri yok (Hale ve Mieher).

**Anlamı:**
- Merkez atom, bağı farklı ve sıfırdan büyük olan bütün komşularını yalnız global alanlarla evrensel olarak yönetebiliyor.
  Tek tek seçilemeyen, birbirine sıkışmış komşular da buna dahil.
- Sınır artık seçicilik değil zaman: süre ~ 1/(bağ farkı). En güçlü 10 komşu için ~ms, yüzlerce komşu için saniyeler mertebesi.
  Bu yüzden elektronun eşevreliliği bu süreyi aşmalı (²⁸Si ya da dinamik ayrıştırma).
- 4'ten büyük kümelerde tam rank sayısal olarak denetlenemedi (boyut 4^(N+1)). Mekanizmanın kendisi her N için sınandı;
  aşağıya bakın.

**Her N için mekanizma (türetme):** ad²_{H₀}(Σ_k I_x^k) = ¼ Σ_k A_k² I_x^k; yani 2j katlı komütatör ¼ʲ Σ_k A_k^{2j} I_x^k verir.
A_k'lar farklı ve sıfırdan büyükse bu katsayılar bir Vandermonde matrisi oluşturur ve matris tersinirdir. O zaman her I_x^k
global komütatörlerin doğrusal birleşimi olur, yani her komşu tek başına sürülebilir. Tek tek çekirdek dönmeleri ve S_z I_z^k
bağlarıyla (yıldız grafiği, bağlı) cebir tam su(2^{N+1}) olur.

Sayısal sınama (`control_all.json`; UN1–UN3 koşudan önce kaydedildi; tam matrisler, B ∥ [1,2,3]):
- N = 6: rank 6/6, çekirdek 1'in I_x'i 7·10⁻¹⁴ hatayla kuruluyor.
- N = 8: rank 8/8, hata 1,6·10⁻¹².
- N = 10: kurma hatası 6·10⁻⁹ (< 10⁻⁶ ✓), ama sayısal rank 8/10 (tekil değer oranı 5·10⁻¹⁷, çift duyarlık sınırı). UN1 burada
  kısmen ✗. Fiziksel anlamı: çok sayıda yakın bağda kontrol aşırı katılaşıyor ve süre hızla uzuyor.
- Denetim (A₁ = A₂): rank 5/6, hata 0,71 ✓ (eşit bağlı iki komşu ayrılamaz).
- Kesin sınama (`vandermonde_exact.json`). Vandermonde determinantı Π_{i<j}(A_j² − A_i²) olduğundan rank, ancak ve ancak bütün
  A_k² farklıysa tamdır; kesirli (tam) aritmetikle denetlendi:
  - N = 10 ve 20'de eşit çift yok, rank kesin olarak tam. N = 10'daki 8/10 yalnızca çift duyarlıklı hesabın sınırıydı.
  - N = 40'ta 1 eşit çift var; ters simetri eşi (r, −r), modelin yapaylığı.
  - N = 100'de 3 eşit çift var; 2'si ters simetri eşi, 1'i genel yönlü alanda tesadüfi eşitlik.

**Gerçekçi eşevrelilikle "bütün komşular":** 6 nm içindeki bütün ²⁹Si'leri ayırmak için gereken süre
T = 1/min(en küçük bağ farkı, en küçük bağ). Ölçülmüş elektron T2'si: ²⁸Si'de 2 s (Tyryshkin vd., Nature Mater. 11, 143, 2012).

| malzeme | 6 nm içindeki ²⁹Si | medyan T | T < 0,2 s | T > 1 s |
|---|---|---|---|---|
| 50 ppm ²⁹Si (Tyryshkin'in malzemesi) | 2,3 | 0,67 ms | %98 | %2 |
| %0,1 ²⁹Si | 45 | 0,35 s | %36 (UN2 ✗, öngörü ≥ %80) | %25 |
| doğal Si (%4,685) | 2119 | ∞ | %0 | %100 (UN3 ✓) |

Doğal Si'de yapılanmaların hepsinde, %0,1'de %3'ünde birebir eşit bağ var; bunlar modelin simetri eşitlikleri.

**Sonuç:**
- Bütün spinli komşuları (6 nm içi) ölçülmüş eşevrelilik süresi içinde evrensel kontrol etmek, ²⁹Si'si seyreltilmiş
  silisyumda mümkün:
  - 50 ppm'de yapılanmaların %98'inde, 0,2 s'den kısa sürede;
  - %0,1'de yapılanmaların %75'inde, < 1 s (T2e'nin yarısı) içinde.
- Doğal silisyumda mümkün değil: yüzlerce komşu, birebir eşit bağlar, saniyeleri aşan süreler.
- Spinsiz ²⁸Si çekirdeklerinin kontrol edilebilecek bir spin serbestlik derecesi yok. Onlarla etkileşim yalnızca elektrostatik
  (alan, §2) ve gerilme yoluyla olur.

## Sonuç

- **Programlanabilir atom evrenselleştirilebiliyor (model düzeyinde):**
  - klasik olarak, alanla bağlı üç atom NAND kurar ve her atomun alanı bütün yakın komşularına ulaşır;
  - kuantum olarak, merkez atom bütün yakın çekirdeklerle etkileşir; ~5 tanesini tek tek, tam kapı setiyle yönetip
    elektronla birlikte 5 parçacıklı dolanıklık kurar.
- **"Bütün komşular":**
  - Etkileşim bütün komşularla var.
  - Tek darbeyle seçerek kontrol doğal Si'de ~5 çekirdekle sınırlı (spektral kalabalık).
  - Ama merkez atom üzerinden global kontrolle, bağı farklı ve sıfırdan büyük olan bütün komşular evrensel olarak
    kontrol edilebilir (§4): Lie cebiri tam (4 çekirdekli 10/10 küme, genel yönlü B); sıkışık komşuya CNOT 0,999–1,000;
    3 çekirdekli GHZ 1,000.
  - Sınır zaman (~1/bağ farkı) ve eşevrelilik.
- **Enerji:** doğrudan atom oynatmakta değil, atomun enerji düzeylerini seçmekte (Stark) ve atom başına maliyeti Landauer'e
  yakın tutmakta işe yarıyor.
- **Ölçüm yok.** Önerilen sınamalar:
  - SKY130'a yakın bir yapıda 4 K'de tek atom yük basamakları ve iki atomlu Coulomb kaymaları;
  - ²⁸Si üzerinde tek donörün ENDOR/NMR ile ²⁹Si kaydı.
