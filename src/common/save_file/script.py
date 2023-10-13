## VIASH START
par = {
  'input': 'A multiline\nstring',
  'output': 'output.txt'
}
## VIASH END

# write input to output
with open(par['output'], 'w') as f:
  f.write(par['input'])
