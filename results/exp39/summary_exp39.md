# Exp39 — Denklem 1 (tam bütün) ile teorik iyon-itkili uzay gemisi

Tarih: 2026-10-19. Bu çalışma yalnız teorik/sayısal modeldir; laboratuvar ölçümü veya fiziksel motor sonucu değildir.

## Girdi ve adlandırma

Kullanıcının `Denklemler` klasöründeki üç görseli ayrı iki denklem değildir; tek bir **Denklem 1 grubudur (tam bütün)**:

- `sqrt(1 - v^2/c^2) = 0,1`
- `E = -mc^2` ve `E = hν, f(r_s)` notları
- `sqrt(1 - v^2/c^2) = m0 c^2/E = m0 c^2/(hν) = 0,1`

Simülasyonda sabit 0,1 değeri genelleştirildi:

`ε_n = sqrt(1 - β^2) = m0 c^2/E = m0 c^2/(hν) = 0,1^n = 10^-n`

Dolayısıyla Lorentz katsayısı `γ_n = 1/ε_n = 10^n` olur. Kısa ölçekte 1 nonilyon `10^30` olduğundan katsayı **n=30** basamağında nonilyona ulaşır. Burada 0,1 hızın kendisi değildir; hız `β = sqrt(1 - ε^2)` ile hesaplanır ve c'yi aşmaz.

Öngörüler simülasyon koşulmadan önce `predictions_exp39.json` içinde D1–D7 olarak kaydedildi; yedisi de PASS verdi.

## Referans gemi ve motor hatları

Exp28/Exp34 çapası korundu: 10 kg ıslak kütle, 8 kg kuru kütle, 2 kg yakıt ve 100 W güç. Exp34'ten alınan egzoz hatları:

| Hat | Egzoz hızı | Isp | F/P |
|---|---:|---:|---:|
| EMI fizibil | 3,809 km/s | 388 s | 392,7 µN/W |
| In diz noktası | 44,818 km/s | 4.570 s | 39,95 µN/W |
| In 5 kV formülü | 91,670 km/s | 9.348 s | 19,61 µN/W |
| Foton alt sınırı | c | 30.570.323 s | 0,00334 µN/W |

Relativistik roket için `R = M_ıslak/M_kuru = exp((c/u)φ)` kullanıldı; burada `φ=acosh(γ)` rapiditedir ve `u` egzoz hızıdır. Büyük n değerlerinde `1-β`, kayan nokta iptalini önlemek için `ε²/(1+β)` biçimiyle hesaplandı.

## Ana sonuçlar

### Nonilyon merdiveni

| n | γ | β | 1−β | 8 kg gemi kinetik enerjisi | Alfa Centauri gemi öz süresi* |
|---:|---:|---:|---:|---:|---:|
| 1 | 10 | 0,9949874371 | 5,013×10^-3 | 6,471×10^18 J | 0,4392 yıl |
| 2 | 100 | 0,9999499987 | 5,000×10^-5 | 7,118×10^19 J | 0,04370 yıl |
| 10 | 10^10 | 1−5×10^-21 | 5,000×10^-21 | 7,190×10^27 J | 4,37×10^-10 yıl |
| 30 | 10^30 = 1 nonilyon | 1−5×10^-61 | 5,000×10^-61 | 7,190×10^47 J | 4,37×10^-30 yıl |

\* Yalnız sabit hızlı seyir varsayımıdır; hızlanma ve yavaşlama süreleri/enerjisi ayrıca gerekir. Dünya çerçevesindeki seyir süresi n=30'da da yaklaşık 4,37 yıldır.

İlk basamakta bile 8 kg kuru geminin kinetik enerjisi `6,471×10^18 J` olur. Kayıpsız 100 W kaynakla yalnız bu enerjiyi toplamak `2,051×10^9 yıl` sürer. n=30'da ideal tedarik süresi `2,278×10^38 yıl` olur.

### Exp34 gemisinin gerçek yanışı

10 kg'dan 8 kg'a relativistik roket yanışı 2.001 kütle adımıyla simüle edildi:

| Hat | Yanış süresi | Son Δv | Son γ−1 | Son ε |
|---|---:|---:|---:|---:|
| EMI fizibil | 0,00615 yıl | 0,850 km/s | 4,02×10^-12 | 0,999999999996 |
| In diz | 0,711 yıl | 10,001 km/s | 5,56×10^-10 | 0,999999999444 |
| In 5 kV | 2,963 yıl | 20,455 km/s | 2,33×10^-9 | 0,999999997672 |

Bu nedenle referans iyon gemisi **n=1 basamağına bile ulaşmaz**. n=1 için gereken kütle oranının `log10R` değeri EMI'de 102.324, In dizde 8.696 ve In 5 kV'de 4.251'dir. Egzozun ışık olduğu kuramsal foton alt sınırında bile oran `R=19,95` gerekir; referans gemide oran yalnız `10/8=1,25`tir.

### Frekans ve Planck sınırı

Tek In iyonu için durgunluk frekansı `ν0 = m_In c²/h = 2,586×10^25 Hz` ve dalga boyu `1,159×10^-17 m`dir. Merdivende `ν=γν0`, `λ=h/(γ m_In c)` olur. Dalga boyu Planck uzunluğuna **n=17,856**'da iner; n=17'de Planck boyunun 7,17 katı, n=18'de 0,717 katıdır. Böylece tek-kuantum `hν` okuması nonilyon basamağından önce standart yarı-klasik alanın dışına çıkar.

### Negatif dal ve `f(r_s)` denetimi

Pozitif dal `E=+γm0c²`, denklem grubundaki `m0c²/E=10^-n` eşitliğini tam sağlar. Buna karşılık yazılı `E_-=-γm0c²` dalı `m0c²/E_-=-10^-n` verir. Pozitif kütle ve pozitif frekansla aynı anda `+0,1^n` eşitliği kurulamaz; negatif dal ancak E'nin büyüklük olarak yeniden tanımlanmasıyla veya negatif frekans/kütle varsayımıyla ayrıca formüle edilebilir. n=1'de negatif kinetik eşdeğeri −72 kg'dır; bu, mevcut gösterimle fiziksel bir itki kazancı olarak yorumlanmadı.

`f(r_s)` görsellerde tanımlı değildir. Yalnız benzerlik olarak Schwarzschild kızıla kayma faktörü `ε_g=sqrt(1-r_s/r)` kullanıldı; n=1 için `r/r_s=1,0101`, n büyüdükçe oran 1'e yaklaşır. Bu satır bir sürücü türetimi değildir.

### Alternatif literal okuma

Eğer amaç γ'nin nonilyona ulaşması değil de **n'nin kendisinin 10^30 yapılması** ise sonuç ancak logaritmik tanımlayıcılarla yazılabilir: `log10γ=10^30`, `log10ε=-10^30`, rapidite yaklaşık `2,303×10^30` ve `log10(1-β)≈-2×10^30`. β bu durumda da 1'den küçük kalır; hız c'yi aşmaz. Simülasyon bunu adım adım yinelemedi.

## Hüküm ve sınırlar

Denklem 1'in pozitif dalı, özel görelilikte tutarlı bir hız/enerji merdiveni tanımlar: katsayı `10^n` büyürken hız c'ye asimptotik yaklaşır ve geminin öz seyir süresi kısalır. Ancak aynı hesap, hızlanmayı geminin kendi iyon yakıtıyla yapmanın kütle-oranı ve enerji bedelini astronomik düzeye çıkarır. Exp34 referans gemisi ilk basamağın çok altında kalır. Negatif enerji dalı mevcut biçimiyle işaret çelişkisi taşır; `f(r_s)` tanımlanmadan itki modeline bağlanamaz.

Bu sonuç bir motorun çalıştığını göstermez. Foton hattı yalnız kuramsal alt sınırdır. Hücre beslemesi, nano-aralık rejimi, hızlanma ivmesi, yıldızlararası ortam/radyasyon etkileri ve laboratuvar itki ölçümü bu deneyin kapsamı dışındadır.

## Dosyalar

- `results/exp39/predictions_exp39.json` — koşudan önce kayıtlı D1–D7 ölçütleri
- `results/exp39/ship39.json` — 0–30 sürekli merdiven, tam sayı basamaklar ve yanış yörüngeleri
- `results/exp39/fig_exp39.png` — dört panelli istatistik figürü (görsel okunabilirlik denetiminden geçti)
- `experiments/exp39_nonillion_ship/ship39.py` — simülasyon
- `experiments/exp39_nonillion_ship/check39.py` — bağımsız cebir denetimi (tüm kontroller True)
- `experiments/exp39_nonillion_ship/fig39.py` — figür üretimi
