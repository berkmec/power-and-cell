# Exp35 — Evrensel hücre–atom kontrol yığını ve tam kontrol sınırları

Tarih: 2026-10-17. Kod: `experiments/exp35_universal_control/`
- `uc35.py`: düzeltilmiş v2 çok atomlu desen yazma + bağımsız-bit ölçekleme koşusu;
- `uc35_v3.py`: küresel çalışma noktası düzeltimi ve tavlama destekli kapalı çevrim;
- `fig35.py`: `fig_exp35.png` istatistik figürü (4500×3000 px).

Öngörüler ve düzeltme kayıtları (hiçbiri sonradan düzenlenmedi):
- `predictions_exp35.json`: ilk UC1–UC6 kaydı;
- `predictions_exp35_v2.json`: MONO yönü ve örnekleyici düzeltmesi; v1 çıktısı `universal_control35_v1_invalid.json` olarak arşivlendi ve hüküm için kullanılmaz;
- `predictions_exp35_v3.json`: v2 tanısı üzerine küresel düzeltim + tavlama takibi.

Sonuçlar:
- `universal_control35.json`: düzeltilmiş v2 (G, L10/L20/L40, CL20, ölçekleme ve skor kartı);
- `universal_control35_v3.json`: v3 (G, L20M, CL20M, CL20A);
- `fig_exp35.png`: desen başarısı, zor desen kalitesi, kayıt ölçekleme ve arka-plan bozulması panelleri.

## Soru ve tanımlar

Bu deney, hücreleri ve atomları tek bir fiziksel düğmeyle değil, ortak bir kontrol protokolüyle evrenselleştirmeyi sınar: **sayım/kalibrasyon → nitelikli hedef seçimi → koruma bantlı yazma → okuma/doğrulama → yerel düzeltme → tutma/dondurma → ancak sonra bağlama**.

Burada “tam kontrol”, sıfır hata iddiası değildir. Seçilmiş hedeflerin tamamının önceden ilan edilen hata bütçesi içinde seçilebilmesi, ayarlanabilmesi, okunabilmesi ve tutulabilmesi; hedef dışı bozulmanın da aynı bütçe altında kalması demektir.

Ana hesap, Exp24 Coulomb-camı modelinde 4 K'de çok atomlu bit deseni yazmaktır. Bit eşlemesi Exp24 ile aynıdır: bit 0 = iyonlaşmış B⁻ (n=1), bit 1 = nötr B⁰ (n=0). Kayıtlı hücre penceresinde en az 3 nitelikli atomu olan yedi 60 nm ada (tohum 0, 1, 2, 5, 6, 7, 8), her adada kalibre adres sırasına göre en fazla 4 atom ve üç desen sınıfı (MONO/ALT/RAND) kullanıldı: toplam 21 görev, 72 görev-biti.

Yerel Gauss alan kipleri idealleştirilmiş dilimli-kapı/Stark alan modlarıdır; SKY130'un atom başına kapı çizebildiği iddia edilmez. Kapalı çevrim, hedef ortalama doluluklarını analog okuma gibi kullanır; tek-atış okuma gürültüsü ve kalibrasyon kayması bu deneyde modellenmemiştir. Fiziksel laboratuvar kontrolü yapılmadı.

## Kısa sonuç

**Tek küresel sürüş keyfî desen yazamaz; küresel düzeltimi de içeren kapalı çevrim, bu modelde seçilmiş 3–4 atomluk bütün desenleri modal olarak doğru yazar.**

- Küresel sürüş G, eşik uyumlu MONO'da 7/7; ALT ve RAND'da 0/7'dir. Bu, tek skaler eşiğin adres sırasına göre yalnız basamak desenleri üretebilmesinin doğrudan sayısal karşılığıdır (UC1 ✓).
- v2'de sabit küresel çalışma noktasına eklenen yerel alanlar ve ilk kapalı çevrim yetmedi: zor desenlerde CL20 marjinali 0,756, tam desen anlık oranı 0,357, modal başarı 5/14 (UC2–UC4 ✗).
- v3'te açık-çevrim adres modeli (L20M) de yetmedi: zor desenlerde modal 2/14. Buna karşılık küresel düzeltimli kapalı çevrim CL20M bütün görevlerde 21/21 modal doğru sonuca ulaştı; ALT+RAND marjinali 0,992 ve tam desen anlık oranı 0,973 oldu.
- Tavlama (CL20A) ek modal kazanç getirmedi: yine 21/21, marjinal 0,994 ve tam desen 0,980. Dolayısıyla kazancı sağlayan ana unsur, adrese göre ayarlanmış başlangıçtan çok **küresel çalışma noktası + yaz-oku-düzelt geri beslemesidir**.
- Bedel arka plandadır: CL20M/CL20A hedef dışı modal çevrilme oranı bütün görevlerde %6,25/%6,29 ile kayıtlı %5 sınırını aşar (UC3.5 ✗). Bu yüzden sonuç, bütün gövdenin tam izolasyonu değil, **seçilmiş hedef altkümesinde tam kontrole ulaşılmasıdır**.

## Protokol sonuçları (21 görev)

| Protokol | Ortalama bit doğruluğu | Tam desen anlık oranı | Modal tam desen | Hedef dışı çevrilme |
|---|---:|---:|---:|---:|
| v2 G (tek küresel V_B) | 0,751 | 0,325 | 7/21 | tanım gereği 0 |
| v2 L10 (yerel, ℓ=10 nm) | 0,813 | 0,523 | 11/21 | %1,35 |
| v2 L20 (yerel, ℓ=20 nm) | 0,789 | 0,449 | 10/21 | %2,36 |
| v2 L40 (yerel, ℓ=40 nm) | 0,745 | 0,308 | 6/21 | %2,16 |
| v2 CL20 (geri besleme) | 0,837 | 0,571 | 12/21 | %2,62 |
| v3 L20M (medyan V_B + adres modeli) | 0,646 | 0,110 | 2/21 | %5,03 |
| v3 CL20M (+ küresel düzeltim) | 0,994 | 0,981 | **21/21** | %6,25 |
| v3 CL20A (+ tavlama) | 0,996 | 0,985 | **21/21** | %6,29 |

Desen kırılımı belirleyicidir: G, MONO'da 7/7 iken ALT ve RAND'da 0/7'dir. CL20M ve CL20A üç sınıfta da 7/7 modal başarı verir. CL20M'de medyan geri besleme güncelleme sayısı 2'dir; son yerel genliklerin medyanı 18,2 mV, %90 dilimi 43,6 mV ve ±60 mV rayına değen oranı %4,17'dir. Son küresel gerilimler 0,904–0,983 V aralığında, hücre penceresi içindedir.

L20M'nin başarısızlığı ayrıca önemlidir: kalibre adreslerden hesaplanan açık-çevrim yerel kaymalar, etkileşimli gövdenin gerçek eşiklerini yeterince tahmin edemez (zor desenlerde modal 2/14). Evrensel yığının kritik katmanı modelin kendisi değil, ölçümle kapatılan çevrimdir.

## Hükümler

### v2 (UC1–UC6)

- **UC1 ✓ (küresel sınır).** MONO 7/7, ALT 0/7, RAND 0/7; G ortalama bit doğruluğu 0,751.
- **UC2 ✗ (yerel kanallar).** ALT+RAND modal başarı G'de 0/14, L10'da 4/14, L20'de 3/14 ve L40'ta 0/14; L20'nin istenen +30 puan üstünlüğü yoktur.
- **UC3 ✗ (ilk kapalı çevrim).** CL20 zor desenlerde marjinal 0,756, tam desen 0,357 ve modal 5/14 verir.
- **UC4 ✗ (korelasyon açığı biçimi).** G'de marjinal çarpımı ile tam desen ortalaması pratikte aynıdır (fark 0,000); başarısızlık çarpım açığından değil, marjinallerin desenle uyumsuz eşik sırasındandır.
- **UC5 ✓ (bağımsız-bit ölçekleme referansı).** Kalibreli havuzda tek atış tam kayıt N=10'da 0,916, N=100'de 0,420; en fazla 3 mükemmel-okumalı denemeyle N=100'de 0,990. Kalibresiz havuzda N=100 tek atış 1,68·10⁻⁹ ve 3 denemede yalnız 0,0047'dir.
- **UC6 ✗ (mesafe/çapraz konuşma biçimi).** L20 yakın-uzak farkı yalnız 0,027 çıktı; CL20 farkı 0,007 ile küçük ama ilk koşul tutmadığı için hüküm faildir.

### v3 (UC3.1–UC3.6)

- **UC3.1 ✗ (açık-çevrim adres modeli).** L20M ALT+RAND modal 2/14 ve bit doğruluğu 0,687.
- **UC3.2 ✓ (küresel düzeltimli geri besleme).** CL20M ALT+RAND'da marjinal 0,992, tam desen 0,973, modal 14/14; medyan 2 güncelleme.
- **UC3.3 ✗ (tavlamanın ek üstünlüğü).** CL20A asıl kalite eşiklerini geçer (0,994; 0,980; 14/14) fakat CL20M de tavana ulaştığı için kayıtlı “en az +10 puan modal üstünlük” koşulu tutmaz.
- **UC3.4 ✓ (ALT yalıtımı).** CL20A, ALT sınıfında 7/7 modal başarı verir (v2 CL20: 0/7).
- **UC3.5 ✗ (arka-plan bedeli).** CL20A hedef dışı çevrilme %6,29 ile %5 sınırının üstündedir; CL20M de %6,25'tir.
- **UC3.6 ✓ (kontrol çabası).** CL20A genliklerinin %90,3'ü ±45 mV içindedir ve bütün son V_B değerleri hücre penceresi içindedir.

## Kayıt ölçekleme dersi

Bağımsız-bit bootstrap yalnızca bir üst/referans hesabıdır; Coulomb-camı koşusu korelasyonları ayrıca gösterir. Yine de kalibrasyonun ölçekteki rolü keskindir:

| Havuz | N=10 tek atış | N=100 tek atış | N=100, ≤3 deneme | N=300, ≤3 deneme |
|---|---:|---:|---:|---:|
| Kalibreli 60 nm (ortalama P=0,9913) | 0,916 | 0,420 | 0,990 | 0,971 |
| Kalibresiz 100 nm, 4 K (ortalama P=0,8157) | 0,130 | 1,68·10⁻⁹ | 0,0047 | 1,05·10⁻⁷ |

Bağımsız lojistik bitlerde bütün kayıt başarısı bit olasılıklarının çarpımıdır; ortak bit marjı için gereken logit marjı yaklaşık log(N/(1−S)) ile büyür. Örneğin %99 tam kayıt için bit başına gereken olasılık N=100'de 0,99990 (logit 9,21), N=1000'de 0,99999'dur (logit 11,51). Geri besleme bu çarpım duvarını ancak okuma hatası, kayma ve etkileşimler izin verdiğinde aşabilir; Exp35'teki mükemmel analog okuma bu yüzden iyimser bir üst sınırdır.

## Katmanlar arası tam-kontrol skor kartı

- **Yük atomları (60 nm gövde):** Adresleme ve tek atom ayarı güçlüdür; çok atomlu keyfî desen bu deneyde yalnız kapalı çevrimle modal olarak tamamlanır. Okuma hâlâ öngörülen akım basamağıdır, ölçülmüş tek-atış dedektörü değildir.
- **İyonik O-boşluğu filamanı:** Kanal sayısı akımla programlanabilir, atom konumları tek tek yerleştirilemez; “atom kontrolü” burada atomik genişlik sayısıdır, konum kaydı değildir.
- **Silisyum filaman hücreleri:** Doğum yeri ve hücre sayısı programlanabilir; tek küresel adımda kazanan kimlikleri güvenilir değildir. Ardışık doğum-bölme-doğrula onarımı gerekir.
- **Saf silisyum spinleri:** Saklı Exp25/26 sonuçlarına göre adresleme + kontrol süresi + saklamayı küçük N'de birlikte geçebilen tek spin rotası saflaştırılmış silisyumdur.
- **Doğal silisyum spinleri:** Lie-cebri/adresleme tek başına yetmez; en yakın komşu ²⁹Si çiftleri saklama sızıntısıyla tam kontrolü bozar.

Bu yüzden elde edilen evrensellik **protokol ve hata-bütçesi düzeyindedir**: normalize adres/marjin koordinatları, nitelikli altküme seçimi, geri besleme ve dondurma bütün katmanlarda ortaktır; fiziksel düğme, okuma ve bağlama katmana göre değişir.

## Yöntem ve köken notları

- Örnekleyici, Exp24'ün 1/2 tek-çevirme + 1/2 dolu→boş sıçrama öneri dağılımı korunarak korunan dolu/boş listeleriyle uygulanmıştır. Tam koşudan önce Exp24 özgün Metropolis'iyle A/B kabul kapısından geçti (6 denetimin 5'inde |ΔP|≤0,25; kalan nokta yavaş camsı geçişte tohum duyarlılığı gösterdi ve ayrı başlatma tanısıyla incelendi).
- Son değerlendirme 1200 tarama, geri besleme ölçümleri 600 taramadır; sıcaklık 4 K'dir. Adresler yeni kalibrasyondan değil, saklı Exp24 denge kalibrasyonlarından alınmıştır.
- `universal_control35.json` içindeki `run_config.note` alanı eski “sparse-hop” metnini korur; v2 hükmü üreten fiilî kod düzeltilmiş Exp24-faithful örnekleyicidir. JSON dosyası köken kaydı olduğu için sonradan değiştirilmemiştir.
- v1 çıktısı yalnız arşivdir; MONO yönü ters ve seyrek sıçrama önerisi yetersiz karıştığı için hiçbir v1 desen hükmü kullanılmamalıdır.

## Açık sınırlar ve sıradaki iş

1. **Arka-plan izolasyonu:** Hedef altkümesi tamamlanırken gövdenin %6'dan fazlası modal olarak çevriliyor. Sonraki denetleyici, hedef doğruluğu ile arka-plan bozulmasını birlikte cezalandıran çok amaçlı maliyet ve hedefler arası etkileşim grafiğine göre zamanlama kullanmalıdır.
2. **Gerçek okuma:** Analog ortalama doluluk yerine tek-atış ikili okuma; okuma hatası, örnekleme süresi ve kalibrasyon kayması açıkça modellenmelidir.
3. **Fiziksel yerel alan:** Gauss kiplerinin dilimli kapı/hücre cepheleriyle üretilebilirliği ve ölçülen çapraz konuşma çekirdeği ayrıca sınanmalıdır.
4. **Tutma ve bağlama:** Desen yazıldıktan sonra sürüş azaltıldığında saklama, NAND/koşullu işlem ve hücreden atoma enerji aktarımı aynı görevde birleştirilmelidir.
5. **Laboratuvar:** Exp33 protokolündeki ileri-VB ön ölçümü ve geniş/balastlı SKY130 koşusu, bu modeldeki kapı/aralık varsayımlarını fiziksel girdilere bağlamak için hâlâ gereklidir.

**Sonuç cümlesi:** Tam kontrol, tek küresel sürüşle değil; kalibrasyon, yerel kanallar, küresel çalışma noktası ve yaz-oku-düzelt çevrimi birleştiğinde bu modelde seçilmiş atom altkümeleri için elde edilmiştir. Bütün gövdenin bozulmadan kontrolü ve fiziksel okuma/yerel alan gerçeklemesi hâlâ açıktır.

## Kaynaklar

Exp18–20 (`summary_exp18.md`, `summary_exp19.md`, `summary_exp20.md`), Exp24 (`summary_exp24.md`, `small_bodies.json`, `thermal.json`, `ionic.json`), Exp25–26 (`summary_exp25.md`, `summary_exp26.md`, `control_all.json`, `register_isotopes.json`, `vmrev_survey.json`), Exp33–34 özetleri ve bu deneyin üç öngörü dosyası.

