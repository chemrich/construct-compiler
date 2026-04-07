import argparse
import random
import yaml
from pathlib import Path
import os
import requests
import json
import re

# Lists of components for high variability
GENES = ["GFP", "RFP", "mCherry", "EGFP", "mScarlet", "Luciferase", "BFP", "YFP", "lacZ", "ampR"]
TAGS = ["6xHis", "Flag", "HA", "Myc", "Strep", "MBP", "GST", "SUMO"]
PROMOTERS = ["T7", "T7lac", "lac", "tac", "trc", "pBAD", "pTet", "pLac", "CMV", "SV40"]
HOSTS = ["e_coli", "e_coli_bl21", "e_coli_k12", "s_cerevisiae", "mammalian", "hek293"]
LINKERS = ["GS_flexible", "rigid", "alpha_helical"]

TEMPLATES = [
    "Express {gene} in {host} with a {tag} tag.",
    "I need a construct with {promoter} driving {gene} fused to {tag}.",
    "Design a plasmid for expressing {gene} in {host} using {promoter} promoter.",
    "Create a fusion of {gene} and {tag} in {host}.",
    "Express {gene} with N-terminal {tag} in {host}.",
    "Express {gene} with C-terminal {tag} in {host}.",
    "Polycistronic construct expressing {gene1} and {gene2} in {host}.",
]

def generate_template_prompt(i):
    """Generate a prompt and its expected outcomes using templates."""
    gene = random.choice(GENES)
    host = random.choice(HOSTS)
    tag = random.choice(TAGS)
    promoter = random.choice(PROMOTERS)
    
    template = random.choice(TEMPLATES)
    
    # Handle templates with multiple genes
    if "{gene1}" in template:
        gene1 = gene
        gene2 = random.choice([g for g in GENES if g != gene])
        prompt = template.format(gene1=gene1, gene2=gene2, host=host, tag=tag, promoter=promoter)
        expect = {
            "host": host,
            "min_cistrons": 2
        }
    else:
        prompt = template.format(gene=gene, host=host, tag=tag, promoter=promoter)
        expect = {
            "host": host,
            "min_cistrons": 1
        }
        
    if "{tag}" in template:
        expect["must_have_parts"] = ["PurificationTag"]
        
    return {
        "id": f"v4_template_{i}",
        "category": "template_generated",
        "difficulty": "medium",
        "prompt": prompt,
        "expect": expect
    }

def generate_llm_prompts_batch(count, model_name):
    """Generate prompts using an LLM (Gemini or Claude)."""
    prompts = []
    chunk_size = 20
    chunks = (count + chunk_size - 1) // chunk_size
    
    for i in range(chunks):
        current_chunk_size = min(chunk_size, count - len(prompts))
        print(f"  Generating chunk {i+1}/{chunks} ({current_chunk_size} prompts)...")
        
        system_prompt = "You are an expert in synthetic biology and genetic engineering. Generate diverse natural language prompts that a researcher might use to specify a genetic construct."
        user_prompt = f"Generate {current_chunk_size} diverse prompts for testing a genetic construct compiler. Each prompt should describe a construct design request. Return the results as a JSON list of objects, where each object has 'prompt' and 'expect' keys. The 'expect' key should be a dict containing expected properties like 'host' (e.g. 'e_coli', 'mammalian'), 'min_cistrons' (integer), and 'must_have_parts' (list of strings like 'PurificationTag', 'CleavageSite', 'SolubilityTag')."
        
        if "gemini" in model_name.lower():
            chunk_prompts = _call_gemini(user_prompt, system_prompt, model_name)
        elif "claude" in model_name.lower():
            chunk_prompts = _call_claude(user_prompt, system_prompt, model_name)
        else:
            raise ValueError(f"Unsupported model: {model_name}")
            
        for p in chunk_prompts:
            p["id"] = f"v4_llm_{len(prompts)}"
            p["category"] = "llm_generated"
            p["difficulty"] = "medium"
            prompts.append(p)
            
        if len(prompts) >= count:
            break
            
    return prompts[:count]

def _call_gemini(prompt, system_prompt, model):
    api_key = os.environ.get("GEMINI_API_KEY") or os.environ.get("GOOGLE_API_KEY")
    if not api_key:
        raise ValueError("Neither GEMINI_API_KEY nor GOOGLE_API_KEY environment variable set")
        
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent?key={api_key}"
    
    full_prompt = f"{system_prompt}\n\n{prompt}"
    
    payload = {
        "contents": [{
            "parts": [{
                "text": full_prompt
            }]
        }],
        "generationConfig": {
            "responseMimeType": "application/json"
        }
    }
    
    try:
        response = requests.post(url, json=payload)
        response.raise_for_status()
        data = response.json()
        
        text = data['candidates'][0]['content']['parts'][0]['text']
        return json.loads(text)
    except (KeyError, IndexError, json.JSONDecodeError, requests.RequestException) as e:
        print(f"Error calling or parsing Gemini: {e}")
        return []

def _call_claude(prompt, system_prompt, model):
    try:
        import anthropic
    except ImportError:
        raise ImportError("anthropic package not installed. Run: pip install anthropic")
        
    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        raise ValueError("ANTHROPIC_API_KEY environment variable not set")
        
    client = anthropic.Anthropic(api_key=api_key)
    
    full_prompt = f"{prompt}\n\nPlease respond with ONLY a valid JSON list of objects."
    
    try:
        response = client.messages.create(
            model=model,
            max_tokens=4096,
            system=system_prompt,
            messages=[{"role": "user", "content": full_prompt}],
        )
        
        text = response.content[0].text
        match = re.search(r'\[.*\]', text, re.DOTALL)
        if match:
            return json.loads(match.group(0))
        return json.loads(text)
    except (json.JSONDecodeError, Exception) as e:
        print(f"Error calling or parsing Claude: {e}")
        return []

def main():
    parser = argparse.ArgumentParser(description="Generate prompt corpus for construct compiler evals.")
    parser.add_argument("--count", type=int, default=1000, help="Number of prompts to generate")
    parser.add_argument("--mode", default="template", choices=["template", "llm"], help="Generation mode")
    parser.add_argument("--model", help="LLM model to use (for llm mode)")
    args = parser.parse_args()
    
    prompts = []
    
    if args.mode == "template":
        print(f"Generating {args.count} prompts using templates...")
        for i in range(args.count):
            prompts.append(generate_template_prompt(i))
    elif args.mode == "llm":
        if not args.model:
            raise ValueError("--model is required for llm mode")
        print(f"Generating {args.count} prompts using LLM ({args.model})...")
        prompts = generate_llm_prompts_batch(args.count, args.model)
        
    output = {"prompts": prompts}
    
    output_path = Path(__file__).parent / "prompt_corpus_v4.yaml"
    with open(output_path, 'w') as f:
        yaml.dump(output, f, sort_keys=False)
        
    print(f"Successfully generated {len(prompts)} prompts to {output_path}")

if __name__ == "__main__":
    main()
