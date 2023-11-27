header = []
content = []
with open('./data/sample_input/public_case_8.in') as f:
    for _ in range(3):
        header.append(next(f))
    for line in f:
        content.append(line)

new_content = []
for c in content:
    if len(c) >= 3:
        c = c.split(" ")
        c[2] = str(int(c[2]) + 200)
        c = ' '.join(c)
        c += '\n'
    new_content.append(c)

with open("./data/sample_input/public_case_8_transformed.in", "w") as file:
    # write to file
    for h in header:
        file.writelines(h)
    for c in new_content:
        file.writelines(c)
